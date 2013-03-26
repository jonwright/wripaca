import matplotlib.pylab as pylab
import geometry
import numpy as np
from ImageD11 import transform, indexing, gv_general
from ImageD11.columnfile import columnfile, colfile_from_hdf
from ImageD11.grain import read_grain_file, write_grain_file, grain
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
        self.s_sq = 0
        self.nobs = 0
        for i in range(3):
            for j in range(3):
                self.gvars.append("ub%d%d"%(i,j))
        #print self.gvars
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
        self.s_sq += (diff*diff).ravel().sum()
        self.nobs += len(diff.ravel())


    def solve(self):
        self.chi2 = self.s_sq / ( self.nobs - self.nv )
        print "Inverting the matrix...", self.chi2
        if 0:
            pylab.clf()
            pylab.imshow( np.log(abs(self.lsq)+1e-10), interpolation='nearest' )
            pylab.show()
        self.imat = np.linalg.inv(self.lsq)
        self.shifts = np.dot( self.imat, self.rhs )
        return self.imat, self.shifts

    def applyshifts( self, pars, grs ):
        errs = np.sqrt(np.diagonal(self.imat)*self.chi2)
        for i,v in enumerate(self.pars):
            print i,v,self.shifts[i],
            pars.parameters[v] += self.shifts[i]
            print "--->",pars.parameters[v],"+/-",errs[i]
        for ig in range(len(grs)):
            offset = self.np + 12*ig
            for i,p in enumerate(self.gvars):
#                print ig,i,p,self.shifts[i+offset]
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
            grs[ig].ubi = np.linalg.inv( grs[ig].ub )
        # max shift over error
        return abs(self.shifts/errs).max()

import wripaca
from ImageD11 import gv_general

def omegacalc_ub( ubi, h, pars, romega):
    pre = np.eye(3).ravel()
    posti = np.dot(gv_general.wedgemat(pars.get('wedge')), 
                     gv_general.chimat(pars.get('chi'))).T.ravel()
    axis = np.array([0,0,-1],np.float)
    ub = np.linalg.inv(ubi)
    gcalc = np.dot( ub, h.T ).T.copy()
    romegacalc = np.zeros( romega.shape, np.float)
    romegaerr  = np.zeros( romega.shape, np.float)
    wripaca.omegacalcclose(gcalc,
                           pre,
                           posti,
                           axis,
                           romega,
                           romegacalc,
                           romegaerr,
                           pars.get('wavelength'),
                           )
    return 180*romegacalc/np.pi, 180*romegaerr/np.pi



def weightsfromdiff( diff, nmedians=2):
    assert diff.shape[0] == 3
    ng = diff.shape[1]
    modlen = np.sqrt( (diff*diff).sum(axis=0) )
    middle = np.median(modlen)
    norm = modlen/(middle*nmedians)
    wt = np.exp( -norm*norm )
    return wt
    if 0:
        pylab.clf()
        pylab.subplot(121)
        pylab.hist( modlen, bins=50 )
        pylab.subplot(122)
        pylab.plot( modlen, wt, ",")
        pylab.show()


def fitgrainfunc( cf, gr, pars, refinables, omegaslop = 0.1 ):
#    print cf.omegaobs[:10]
    omc, ome = omegacalc_ub( gr.ubi, cf.hkl, pars,
                cf.omegaobs*np.pi/180.0 )
    
    #
    # wt =1 - np.exp(-ome*ome/omegaslop)
    # pylab.plot( ome, wt, ",")
    # pylab.show()
    cf.omega = np.where( ome < omegaslop, omc, cf.omegaobs )
    
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
    

def fitmanygrains( cfs, grs, pars, refinables, symmetrise ):
    mlsq = many_grain_lsq(len(grs), refinables)
    e = 100
    h = 0
    totalerr = 0
    for i, cf, gr in zip(range(len(cfs)), cfs, grs):
        pars.set('t_x', gr.translation[0])
        pars.set('t_y', gr.translation[1])
        pars.set('t_z', gr.translation[2])
        diff, Ddiff = fitgrainfunc( cf, gr, pars, refinables )
        wt = weightsfromdiff( diff )
        np.multiply( diff, wt, diff)
        for p in Ddiff.keys():
#            print p,wt.shape, Ddiff[p].shape, diff.shape
            np.multiply( Ddiff[p], wt, Ddiff[p] )
        drlv = (diff*diff).sum(axis=0)
        hnew,e = np.histogram( drlv, bins=e )
        h = hnew + h
        totalerr += drlv.sum()
        mlsq.addgrain( i, diff, Ddiff )
    imat,shifts = mlsq.solve()
    max_shift = mlsq.applyshifts( pars, grs )
    if symmetrise:
        grs = [ grain( makehexagonal(g.ubi), g.translation) for g in grs ]

    #pylab.plot( (e[1:]+e[:-1])*0.5, np.log(h+0.5), "o-")
    print "totalerr",totalerr, max_shift
    return pars, grs, max_shift



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
    #print "entry to assignpeaks"
    #print "omega[:10]",colfile.omega[:10]
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
        pylab.figure()
        pylab.plot(r_omega - romegacalc, romegaerr,",")
        pylab.show()
    r_omega = np.where( abs(romegaerr) < omegatol*np.pi/180.,  romegacalc, r_omega )
    cyf.hkl( colfile.sc, colfile.fc, r_omega, gr, hkl, kcalc )
    drlv_omegafree = drlv.copy()
    wripaca.ih_drlv( hkl, drlv_omegafree )
    m = drlv_omegafree < tol
    sc = colfile.sc[m]
    fc = colfile.fc[m]
    omegaobs = colfile.omega[m].astype(np.float)
    omega = colfile.omega[m].astype(np.float)
    hkl = hkl[m]
    print "Got",m.sum(),"peaks",
#    print "end of peak selection"
#    print "colfile.omega[:10]",colfile.omega[:10]
#    print "omega[:10]", omega[:10]
#    print "omegaobs[:10]",omegaobs[:10]
    if 0:
        pylab.ion()
        pylab.clf()
        pylab.title("tol = %f"%(tol))
        pylab.hist( drlv, bins=np.arange(0,0.05,0.001))
        pylab.hist( drlv_omegafree, bins=np.arange(0,0.05,0.001))
        pylab.show()
        raw_input("OK?")
    omegaobs = (180./np.pi)*np.array([wripaca.anglediff(o,0.0) for o in
        omegaobs*np.pi/180.])
    return dummycf( sc, fc, omega, omegaobs, hkl )

def fitallgrains( gfile, pfile, cfile, ngrains = None, tolerance=0.02,
        symmetrise = False):
    colfile = loadcolfile( cfile )
    if symmetrise:
        grains = [ grain( makehexagonal(g.ubi), g.translation) for g in
            read_grain_file( gfile ) ]
    else:
        grains = read_grain_file( gfile ) 
    pars = read_par_file( pfile )

    variables = [ 't_x','t_y', 't_z', 'y_center',  'tilt_y', 'tilt_z',
                       'tilt_x',  'distance', 'wedge']
 #   variables = [ 't_x','t_y', 't_z']
    #, 'y_center',  'tilt_y', 'tilt_z',
    #                   'tilt_x',  'distance', 'wedge']
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
        cf = assignpeaks( gr, pars, colfile, tol = tolerance )
        cfs.append(cf)
        grs.append(gr)

    pi = parameters( **pars.parameters )
    refpars, grs, shoemax = fitmanygrains( cfs, grs, pi, variables, symmetrise )
    for i in range(10):
        refpars, grs, shoemax = fitmanygrains( cfs, grs, refpars, variables, symmetrise )
        if shoemax < 1:
            break
    pylab.show()
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


    print "Writing grains to",gfile
    write_grain_file( gfile, grs)
    print "Writing pars to",pfile
    refpars.saveparameters( pfile )
    print indexing.ubitocellpars( grs[0].ubi )


if __name__=="__main__":
    fitallgrains( sys.argv[1], sys.argv[2], sys.argv[3], ngrains = None,
            tolerance=float(sys.argv[4]) )


