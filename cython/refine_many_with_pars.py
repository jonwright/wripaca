
import matplotlib.pylab as pylab
import geometry
import numpy as np
from ImageD11 import transform, indexing, gv_general
from ImageD11.columnfile import columnfile, colfile_from_hdf
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
            pars.parameters[p] += shifts[i]
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

class many_grain_lsq:
    def __init__(self, ngrains, pars):
        self.pars = [p for p in pars] # globals
        self.pars.remove('t_x')
        self.pars.remove('t_y')
        self.pars.remove('t_z')
        self.np = len(self.pars)
        self.nv = 12*ngrains + len(self.pars)
        self.rhs = np.zeros(self.nv)
        self.lsq = np.zeros((self.nv,self.nv))
        self.gvars = ['t_x','t_y','t_z']
        for i in range(3):
            for j in range(3):
                self.gvars.append("ub%d%d"%(i,j))
        print "The size of the least squares problem is", self.nv
        # FIXME this should be block banded

    def addgrain( self, index, diff, Ddiff ):
        # par - par terms
        for i in range(len(self.pars)):
            gi = Ddiff[self.pars[i]]
            self.rhs[i] += (gi*diff).ravel().sum()
            self.lsq[i,i] += (gi*gi).ravel().sum()
            for j in range(i+1,len(self.pars)):
                gj = Ddiff[self.pars[j]]
                t = (gi*gj).ravel().sum()
                self.lsq[i,j] += t
                self.lsq[j,i] += t
        # grain - grain : overwrite
        goffset = self.np + 12*index
        for i in range(len(self.gvars)):
            gi = Ddiff[self.gvars[i]]
            ip = i + goffset 
            self.rhs[ ip ] = (gi*diff).ravel().sum()
            self.lsq[ ip, ip ] = (gi*gi).ravel().sum()
            for j in range(i+1, len(self.gvars)):
                gj = Ddiff[self.gvars[j]]
                jp = j + goffset
                t = (gi*gj).ravel().sum()
                self.lsq[ip,jp] = t
                self.lsq[jp,ip] = t
        # grain - par : overwrite
        goffset = self.np + 12*index
        for i in range(len(self.gvars)):
            gi = Ddiff[self.gvars[i]]
            ip = i + goffset 
            self.rhs[ ip ] = (gi*diff).ravel().sum()
            self.lsq[ ip, ip ] = (gi*gi).ravel().sum()
            for j in range(len(self.pars)):
                gj = Ddiff[self.pars[j]]
                t = (gi*gj).ravel().sum()
                self.lsq[ip,j] = t
                self.lsq[j,ip] = t


    def solve(self):
        print "Inverting the matrix..."
        self.imat = np.linalg.inv(self.lsq)
        self.shifts = np.dot( self.imat, self.rhs )
        return self.imat, self.shifts

    def applyshifts( self, pars, grs ):
        for i,v in enumerate(self.pars):
            print i,v,self.shifts[i],
            pars.parameters[v] += self.shifts[i]
            print "--->",pars.parameters[v]
        for ig in range(len(grs)):
            offset = self.np + 12*ig
            for i,p in enumerate(self.gvars):
                if p[0:2] == "ub":
                    ip = int(p[-2])
                    jp = int(p[-1])
                    grs[ig].ub[ip,jp] += self.shifts[i+offset]
                if p == 't_x':
                    grs[ig].translation[0] += self.shifts[i+offset]
                if p == 't_y':
                    grs[ig].translation[1] += self.shifts[i+offset]
                if p == 't_z':
                    grs[ig].translation[2] += self.shifts[i+offset]








def fitgrainfunc( cf, gr, pars, refinables ):
    cf.omega = (180.0/np.pi)*geometry.omegacalc_ub( gr.ubi, cf.hkl, pars,
            cf.omegaobs*np.pi/180.0 )
    gobs, Dg = geometry.Derivativegobs( cf, gr, pars, refinables )
    gcalc, dGcalc = geometry.derivativegcalc( gr, gobs )
    diff = gobs - gcalc
    Ddiff = dGcalc
    for p in Dg.keys():
        if Ddiff.has_key(p):
            print "wtf"
        else:
            Ddiff[p] = -Dg[p]
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
       # if not quiet:
        print "chi2",chi2,
        shoemax  =applyshifts( shifts, imat, chi2, pars, gr, vars, quiet=quiet )

        if not quiet:
            print "Max shift/e % 6.2g"%(shoemax)

        if 0:
            gr.ubi = makehexagonal( gr.ubi )
        # make W from lsq

        cellparerrors(gr, vars, imat, chi2, quiet=quiet)
        cycle += 1
    return gr, pars


def fitmanygrains( cfs, grs, pars, refinables ):
    mlsq = many_grain_lsq(len(grs), refinables)
    for i, cf, gr in zip(range(len(cfs)), cfs, grs):
        pars.set('t_x', gr.translation[0])
        pars.set('t_y', gr.translation[1])
        pars.set('t_z', gr.translation[2])
        diff, Ddiff = fitgrainfunc( cf, gr, pars, refinables )
        print i, (diff*diff).sum()
        mlsq.addgrain( i, diff, Ddiff )
    imat,shifts = mlsq.solve()
    mlsq.applyshifts( pars, grs )
    return pars

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




# adding labels...
# foreach grain compute drlv with ideal omega
# choose lowest


def loadcolfile(cfile):
    """ fixme: move to ImageD11.colfile if not there already """
    if cfile.find("::")>-1:
        h,s = cfile.split("::")
        try:
            colfile = colfile_from_hdf( h, name = s )
        except:
            print "*",cfile,"*"
            print "Trying to open",h, os.path.exists(h)
            print s
            raise
    else:
        colfile = columnfile(cfile)
    return colfile


import cyfit, wripaca

class dummycf:
    def __init__(self, sc, fc, omega, omegaobs, hkl):
        self.nrows= len(sc)
        self.sc = sc
        assert len(sc) == self.nrows
        self.fc = fc
        assert len(fc) == self.nrows
        self.omega = omega
        assert len(omega) == self.nrows
        self.omegaobs = omegaobs
        assert len(omegaobs) == self.nrows
        self.hkl = hkl
        assert len(hkl) == self.nrows

def assignpeaks( gr, pars, colfile, tol=None, omegatol=0.2 ):
    labels = np.zeros( colfile.nrows )
    cyf = cyfit.cyfit(pars)
    cyf.setscfc( colfile.sc, colfile.fc )
    hkl = np.zeros( cyf.XL.shape, np.float)
    drlv = np.zeros( colfile.nrows, np.float)
    kcalc = np.zeros( cyf.XL.shape, np.float)
    r_omega  = np.array(colfile.omega*np.pi/180.0,np.float)
    cyf.hkl( colfile.sc, colfile.fc, r_omega, gr, hkl, kcalc )
    # Make integer hkl
    wripaca.ih_drlv( hkl, drlv )
    # Now replace omegaobs by omegacalc where this is appropriate
    pre = np.eye(3).ravel()
    posti = np.dot(gv_general.wedgemat(pars.get('wedge')), 
                     gv_general.chimat(pars.get('chi'))).T.ravel()
    axis = np.array([0,0,-1],np.float)
    ub = np.linalg.inv(gr.ubi)
    gcalc = np.dot( ub, hkl.T ).T.copy()
    #print hkl
    romegacalc = np.zeros( r_omega.shape, np.float)
    romegaerr  = np.zeros( r_omega.shape, np.float)
    wripaca.omegacalcclose(gcalc,
                           pre,
                           posti,
                           axis,
                           r_omega,
                           romegacalc,
                           romegaerr,
                           pars.get('wavelength'),
                           )
    # OK, now we will accept anything within omegatol
    if 0:
        pylab.hist(romegaerr, bins=50)
        pylab.figure()
        pylab.plot(r_omega, romegaerr,",")
        pylab.show()
    r_omega = np.where( romegaerr < omegatol*np.pi/180.,  romegacalc, r_omega )
    cyf.hkl( colfile.sc, colfile.fc, r_omega, gr, hkl, kcalc )
    drlv_omegafree = drlv.copy()
    wripaca.ih_drlv( hkl, drlv_omegafree )
    m = drlv < tol
    sc = colfile.sc[m]
    fc = colfile.fc[m]
    omegaobs = colfile.omega[m].astype(np.float)
    omega = colfile.omega[m].astype(np.float)
    hkl = hkl[m]
    print "Got",m.sum(),"peaks",
    if 0:
        pylab.ion()
        pylab.title("tol = %f"%(tol))
        pylab.hist( drlv, bins=np.arange(0,0.05,0.001))
        pylab.hist( drlv_omegafree, bins=np.arange(0,0.05,0.001))
        pylab.show()
        raw_input("OK?")

    return dummycf( sc, fc, omega, omegaobs, hkl )

def fitallgrains( gfile, pfile, cfile, ngrains = None):
    colfile = loadcolfile( cfile )
    grains = read_grain_file( gfile )
    pars = read_par_file( pfile )

    variables = [ 't_x','t_y', 't_z', 'y_center',  'tilt_y', 'tilt_z',
                       'tilt_x',  'distance', 'wedge']
    pfitted = []
    grs = []
    cfs = []
    if ngrains is None:
        ng = len(grains)
    else:
        ng = ngrains
    
    for i in range(ng):
        print "***",i,
        gr = grains[i]
        cf = assignpeaks( gr, pars, colfile, tol = 0.02 )
        cfs.append(cf)
        grs.append(gr)
    pi = parameters( **pars.parameters )
    refpars = fitmanygrains( cfs, grs, pi, variables )
    for i in range(3):
        refpars = fitmanygrains( cfs, grs, refpars, variables )

    if 0:
        pi = parameters( **pars.parameters )
        pi.set('t_x', gr.translation[0])
        pi.set('t_y', gr.translation[1])
        pi.set('t_z', gr.translation[2])
        diff, Ddiff = fitgrainfunc( cf, gr, pi, variables )
        print "%.5g"%((diff*diff).ravel().sum()),
        gr, pfit = fitgrain( cf, gr, pi, variables, quiet=True )
        grains[i] = gr
        pfitted.append( pfit )
        diff, Ddiff = fitgrainfunc( cf, gr, pfit, variables )
        print "%.5g"%((diff*diff).ravel().sum())

        if 0:
            v = Ddiff.keys()
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
    fitallgrains( sys.argv[1], sys.argv[2], sys.argv[3], ngrains = None )


