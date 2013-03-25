

import geometry
import numpy as np
from ImageD11 import transform, indexing
from ImageD11.columnfile import columnfile
from ImageD11.grain import read_grain_file, write_grain_file
from ImageD11.parameters import parameters, read_par_file
from ImageD11.transformer import transformer
from ImageD11.unitcell import unitcell
import sys, time

def makehexagonal( ubi ):
    """ Finds a nearby hexagonal cell """
    cp = indexing.ubitocellpars( ubi )
    a = (cp[0]+cp[1])/2
    cell = unitcell( (a,a,cp[2],90,90,120),"P")
    u = indexing.ubitoU( ubi )
    ub = np.dot( u, cell.B )
    return np.linalg.inv( ub )

def cpfromub( ubf ):
    ub = np.reshape( ubf, (3,3) )
    return np.array(indexing.ubitocellpars( np.linalg.inv( ub ) ))

def applyshifts( shifts, imat, chi2, pars, gr, vars, quiet=True):
    # Apply the shifts
    shoemax = 0.0 # max shift over error
    for i in range(len(vars)):
        p = vars[i]
        e = np.sqrt( imat[i,i]*chi2 )
        if pars.parameters.has_key(p):
            if not quiet:
                print "%-20s %20.9f"%(p, pars.get(p)+shifts[i])
                print "%-20s %20.9f"%(p+"_error", e )
                print "%-20s %20.9f"%(p+"_sh/e", shifts[i]/e )
            pars.parameters[p] -= shifts[i]
        elif p.find("ub")==0:
            ip = int(p[-2])
            jp = int(p[-1])
            gr.ub[ip,jp] += shifts[i]
            if not quiet:
                print "%-20s %20.9f"%(p, gr.ub[ip,jp])
                print "%-20s %20.9f"%(p+"_error", e )
                print "%-20s %20.9f"%(p+"_sh/e", shifts[i]/e )        
        else:
            print "What",p
        shoe = abs(shifts[i] / e)
        if shoe > shoemax:
            shoemax = shoe
    gr.translation = pars.get('t_x'),pars.get('t_y'),pars.get('t_z')
    gr.ubi = np.linalg.inv( gr.ub )
    return shoemax

def cellparerrors(gr, vars, imat, chi2, quiet=True):
    cp = cpfromub( gr.ub.ravel() )
    W = np.zeros((9,9)) + 42
    for i in range(len(vars)):
        pi = vars[i]
        if pi.find("ub")!=0:
            continue
        ni = int(pi[-2])*3 + int(pi[-1])
        for j in range(len(vars)):
            pj = vars[j]
            if pj.find("ub")!=0:
                continue
            nj = int(pj[-2])*3 + int(pj[-1])
            W[ni,nj] = imat[i,j]*chi2
    dcpdub = []
    rv = np.linalg.det( gr.ub ) # recip vol
    #print rv
    ubf = gr.ub.ravel()
    for i in range(9):
        ubfi = gr.ub.copy().ravel()
        s = rv*0.001
        ubfi[i] += s
        #    print (ubfi[i] - ubf[i])/ubf[i]
        cpi = cpfromub( ubfi )
        dcpdub.append((cpi - cp)/s)
    Q = np.array( dcpdub )
    CPE = np.dot( Q.T, np.dot( W, Q ) )
    if not quiet:
        for i in range(6):
            print cp[i],np.sqrt(CPE[i,i])
        print "a==b?:",0.5*(cp[0]-cp[1])/(cp[0]+cp[1])
        print "alpha %e %e"%(np.cos(np.radians(cp[3])), np.radians(cp[3]-90))
        print "beta %e %e"%(np.cos(np.radians(cp[4])),np.radians(cp[4]-90))
        print "gamma %e %e"%(np.cos(np.radians(cp[5]))+0.5,np.radians(cp[5]-120))


def lsq( diff, Ddiff, Weights, allref):
    nv  = len(allref)
    rhs = np.zeros( nv )
    lsq = np.zeros( (nv, nv ) )
    for i in range(nv):
        gi = Ddiff[allref[i]]
        rhs[i] = (gi*diff).ravel().sum() 
        for j in range(i,nv):
            gj = Ddiff[allref[j]]
            t = (gi*gj).ravel().sum()
            lsq[i,j] = t
            lsq[j,i] = t
    imat = np.linalg.inv(lsq)
    shifts = np.dot( imat, rhs )
    chi2 = (diff*diff).ravel().sum()
    return imat, shifts, chi2



def fitgrainfunc( cf, gr, pars, refinables ):
    gobs, Dg = geometry.Derivativegobs( cf, gr, pars, refinables )
    gcalc, dGcalc = geometry.derivativegcalc( gr, gobs )
    diff = gobs - gcalc
    Ddiff = Dg
    for p in dGcalc.keys():
        if Ddiff.has_key(p):
            print "wtf"
        else:
            Ddiff[p] = dGcalc[p]
    return diff, Ddiff
    

def fitgrain( cf, gr, pars, refinables, ncycle=20, shoelim=1e-4, quiet = False):
    """
    Fits a grain to data in colfile using pars and varying variables
    """
    cycle = 0
    shoemax = None
    while cycle < ncycle: # Break when converged
        start = time.clock()
        if shoemax is not None:
            if shoemax < shoelim:
                break
        # FIXME: Omega column
        diff, Ddiff = fitgrainfunc( cf, gr, pars, refinables )
        vars = Ddiff.keys()

        imat, shifts, chi2 = lsq( diff, Ddiff, None, vars )
        print "chi2",chi2,
        shoemax  =applyshifts( shifts, imat, chi2, pars, gr, vars )

        if not quiet:
            print "Max shift/e % 6.2g"%(shoemax)

        if 0:
            gr.ubi = makehexagonal( gr.ubi )
        # make W from lsq

        cellparerrors(gr, vars, imat, chi2, quiet=False)
        cycle += 1
    return gr, pars






# Unit cell esds

# brain dead method





# Now would like to see errors projected onto pixel coords or omega/tth/eta etc

def project_diff_on_variable(diff, variable, Derivs):
    g = Derivs[variable]
    modg = np.sqrt( (g*g).sum(axis=0) )
    ghat = g/modg
    # next is a scalar quantity
    sdiff_along_ghat =  (ghat * diff).sum(axis=0)
    # This difference in units of g
    return sdiff_along_ghat/modg


import matplotlib.pylab as pylab


# adding labels...
# foreach grain compute drlv with ideal omega
# choose lowest

def computelabels( grains, pars, colfile):
    for g in grains:
        pass

def fitallgrains( gfile, pfile, cfile):
    colfile = columnfile(cfile)
    grains = read_grain_file( gfile )
    pars = read_par_file( pfile )
    variables = [ 't_x','t_y', 't_z', 'y_center',  'tilt_y', 'tilt_z',
                       'tilt_x',  'distance', 'wedge']
    pfitted = []
    for i in range(len(grains)):
        print "***",i
        gr = grains[i]
        cf = colfile.copy()
        cf.filter( cf.labels == i )
        print cf.nrows
        pi = parameters( **pars.parameters )
        pi.set('t_x', gr.translation[0])
        pi.set('t_y', gr.translation[1])
        pi.set('t_z', gr.translation[2])
        gr, pfit = fitgrain( cf, gr, pi, variables )
        grains[i] = gr
        pfitted.append( pfit )

        diff, Ddiff = fitgrainfunc( cf, gr, pfit, variables )
        v = Ddiff.keys()
        pd =   project_diff_on_variable( diff, 'y_center', Ddiff)
        print pd.shape, colfile.omega.shape
        if 0:
            for v in ['y_center', 'distance']:
                pylab.figure(1)
                pylab.title("Versus omega")
                pylab.plot( cf.omega, project_diff_on_variable( diff, v, Ddiff) , ",", label=v) 
                pylab.figure(2)
                pylab.title("Versus fc")
                pylab.plot( cf.fc, project_diff_on_variable( diff, v, Ddiff) , ",", label=v) 
                pylab.figure(3)
                pylab.title("Versus sc")
                pylab.plot( cf.sc, project_diff_on_variable( diff, v, Ddiff) , ",", label=v) 
                pylab.figure(4)
                pylab.title("Versus sigo")
                pylab.plot( cf.sigo, project_diff_on_variable( diff, v, Ddiff) , ",", label=v)    
                pylab.legend()
                pylab.show()
                raw_input()



    write_grain_file( gfile+".fit", grains)
    


if __name__=="__main__":
    fitallgrains( sys.argv[1], sys.argv[2], sys.argv[3] )


