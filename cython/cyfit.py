

import sys, time

try: 
    import wripaca
    print wripaca.__file__
except:
    sys.path.append("build/lib.win-amd64-2.7")
    import wripaca
    print wripaca.__file__

from ImageD11.grain import read_grain_file, grain
from ImageD11.parameters import read_par_file
from ImageD11 import transform, gv_general, indexing
import geometry
import h5py, sys, numpy as np
import scipy.optimize

class cyfit(object):
    def __init__(self, p):
        self.pre = np.eye(3).ravel()
        self.pars = p
        self.reset()

    def reset(self):
        self.post = gv_general.chiwedge( wedge = self.pars.get('wedge'),
                chi = self.pars.get('chi')).ravel()
        self.posti= np.dot(gv_general.wedgemat(self.pars.get('wedge')), 
                           gv_general.chimat(self.pars.get('chi'))).T.ravel()

        self.axis = np.array([0,0,-1],np.float)

    def setscfc(self, sc, fc):
        self.sc = sc
        self.fc = fc
        self.XL = transform.compute_xyz_lab( [sc, fc], 
                **self.pars.parameters ).T.copy()

    def hkl(self, sc, fc , r_omega, g, hkl, kcalc):
        wripaca.hklcalc_many( self.XL,
                self.axis,
                r_omega,
                g.ubi.ravel(),
                g.translation,
                self.pre,
                self.post,
                hkl,
                self.pars.get('wavelength'),
                kcalc)
        return hkl




class refiner(object):
    def __init__(self, g, sc, fc, romega, h, kcalc, cyf):
        self.g = g 
        self.sc = sc
        self.fc = fc
        self.romega = romega
        self.romegacalc = np.zeros( romega.shape, np.float)
        self.romegaerr = np.zeros( romega.shape, np.float)
        self.h = h
        self.kcalc = kcalc
        self.cyf = cyf
        self.wt = np.ones(  romega.shape, np.float )
        self.kcalc = np.zeros( (romega.shape[0],3), np.float)
        self.p0 = np.array(list(self.g.ubi.ravel()) + list( self.g.translation ))
        self.diffs = np.zeros( (romega.shape[0],4), np.float)

    def omegacalc_ub(self, p):
        ubi = p[:9].reshape(3,3)
        ub = np.linalg.inv(ubi)
        #print "ub",ub
        #print "ubi",ubi
        self.gcalc = np.dot( ub, self.h.T ).T.copy()
        #print 'gcalc[0]',self.gcalc[0]
        #print 'h[0]',self.h[0]
        wripaca.omegacalcclose(self.gcalc,
                               self.cyf.pre,
                               self.cyf.posti,
                               self.cyf.axis,
                               self.romega,
                               self.romegacalc,
                               self.romegaerr,
                               self.cyf.pars.get('wavelength'),
                               )
        # e == 0 -> sigma = step
        if 0:
            halfstep = 0.1*np.pi/180.
            abserr = abs(self.romegaerr)
            self.wt = np.where( abserr < 2*halfstep ,
                1.0/(abs( halfstep - abserr)+halfstep*0.3),
                0.1/abserr )
        #print 'self.romegaerr[0]',self.romegaerr[0]
        return self.romegaerr * self.wt

    def fit_ubi_omega(self):
        
        print indexing.ubitocellpars(self.g.ubi)
        f = self.omegacalc_ub
        res = scipy.optimize.leastsq( f, self.p0, full_output=1)
        pfit, pcov, info, errmsg, ier = res
        ub = pfit.reshape(3,3)
        ubi = np.linalg.inv(ub)
        print indexing.ubitocellpars(ubi)

    def fit_trans(self, fp):
        ubi = fp[:9].reshape(3,3)
        t = fp[9:12]
        gr = grain( ubi, translation = t )
        self.omegacalc_ub(fp)
        #print "self.romega[0]", self.romega[0]
        wripaca.rotate_vectors_axis_angle(
                self.cyf.axis,
                -self.romega,
                self.gcalc,
                self.kcalc)
        #print self.kcalc[0]
        #print self.kcalc.shape
        #print self.kcalc[0]
        #  np.dot( self.cyf.posti.reshape(3,3), self.kcalc.T )
        #print self.kcalc[0]
        w = self.cyf.pars.get('wavelength')
        np.add(self.kcalc[:,0],1./w,self.kcalc[:,0])
        mk = np.sqrt((self.kcalc*self.kcalc).sum(axis=1))
        np.divide(self.kcalc.T, mk, self.kcalc.T)
        lv=geometry.labvec( self.sc, self.fc, self.romegacalc*180/np.pi,
                gr, self.cyf.pars)
        modlv = np.sqrt((lv*lv).sum(axis=0))
        kobs = lv/modlv
        #print "kobs",kobs.T[0]
        if 0:
            gr.translation[2]=0
            lvt = geometry.labvec( self.sc, self.fc, self.romegacalc*180/np.pi,
                gr, self.cyf.pars)
            modlvt = np.sqrt((lvt*lvt).sum(axis=0))
            kobst = lvt/modlv
            print kobs[:,:3] - kobst[:,:3]
            sys.exit()
        #print "obs",kobs.shape,"calc",kcalc.shape,self.diffs[:,:3].shape
        self.diffs[:,:3] = kobs.T - self.kcalc
        self.diffs[:,3] = self.romegaerr
        #print "diffs",self.diffs[0]
        return self.diffs.ravel().copy()

    def do_fit_trans(self):
        f = self.fit_trans
        p0 = self.p0.copy()
        #        p0[:-3]=0.0
        #print p0,"call d0"
        #d0 = f(p0)
        if 0:
         for i in range(len(p0)):
            pt = p0.copy()
            pt[i] = p0[i]+0.001
            from matplotlib.pylab import clf, ion, title, plot, show
            print pt - p0, pt
            ion()
            clf()
            title("%d"%(i))
            plot(d0, f(pt) - d0, ",")
            show()
            if raw_input()[0] != " ": 
                break


        res = scipy.optimize.leastsq( f, p0, full_output=1)
        pfit, pcov, info, errmsg, ier = res
        if ier not in [1,2,3,4]:
            print s_sq, ier, errmsg
        else:
            residu = f(pfit) 
            s_sq = (residu**2).sum()/(len(residu)-len(p0))
        ubi = pfit[:9].reshape(3,3)
        print ("%.6f "*6)%(indexing.ubitocellpars(ubi))
        print pfit[9:12]
        self.g = grain( ubi, pfit[9:12].copy())

if __name__=="__main__":
    h = h5py.File(sys.argv[1],"r")
    assert sys.argv[2] in h.keys(), h.keys()
    c = h[sys.argv[2]]
    p = read_par_file(sys.argv[3])
    gf = read_grain_file(sys.argv[4])
    tol = float(sys.argv[5])

    sc = c['sc'][:]
    fc = c['fc'][:]
    omega = c['omega'][:]

    cyf = cyfit( p )
    cyf.setscfc( sc, fc )
    hkl = np.zeros( cyf.XL.shape, np.float)
    drlv = np.zeros( cyf.XL.shape[0], np.float)
    kcalc = np.zeros( cyf.XL.shape, np.float)
    rs = []
    r_omega  = np.array(omega*np.pi/180.0,np.float)
    start = time.clock()
    for i in range(5):#len(gf)):
        cyf.hkl( sc, fc, r_omega, gf[i], hkl, kcalc )
        wripaca.ih_drlv( hkl, drlv )
        m = drlv < tol
        rs.append( refiner( gf[i], sc[m], fc[m], r_omega[m], 
            hkl[m], kcalc[m], cyf) )
    if 0:
        from matplotlib.pylab import *
        hist(drlv,bins=np.arange(0,0.5,0.01))
        show()



    print time.clock()-start
    start = time.clock()
    for r in rs:
        r.do_fit_trans()
    print "fitting time",time.clock()-start
    if 1:
        from matplotlib.pylab import *
        r = rs[0]
        print r.g.translation
        plot((180/np.pi)*r.romega, (180/np.pi)* r.romegaerr, ",")
#        plot(r.romegaerr, wd, ",")
        show()
    


    


