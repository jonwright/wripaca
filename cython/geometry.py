
import ImageD11.transform
import numpy as np

detpnames = [  'distance', 
                'y_center', 'z_center',
                'tilt_x', 'tilt_y', 'tilt_z',
                'y_size', 'z_size', 
                 'o11','o12','o21','o22'
                 ]
detpsteps = {
            'distance' : lambda distance: 0.001 * distance,
            'y_center' : lambda centre: 0.01,
            'z_center' : lambda centre: 0.01,
            'tilt_x'   : lambda tilt: 0.01*np.pi/180.0,
            'tilt_y'   : lambda tilt: 0.01*np.pi/180.0,
            'tilt_z'   : lambda tilt: 0.01*np.pi/180.0,
            'y_size'   : lambda size: 1.001*size,
            'z_size'   : lambda size: 1.001*size 
            }

def computeXLYLZL( colfile, pars):
    pks = [ colfile.sc, colfile.fc ]
    dargs = {}
    for p in detpnames:
        dargs[p] = pars.get(p)
    xlylzl = ImageD11.transform.compute_xyz_lab( pks, **dargs )
    return xlylzl

def derivativeXLYLZL( colfile, pars, refinables ):
    pks = [ colfile.sc, colfile.fc ]
    dargs = {}
    for p in detpnames:
        dargs[p] = pars.get(p)
    calc = computeXLYLZL( colfile, pars )
    # parameter step sizes
    derivs = {}
    for p in refinables:
        step = detpsteps[p](dargs[p])
        psave = dargs[p]
        dargs[p] += step
        xlylzl1 = ImageD11.transform.compute_xyz_lab( pks, **dargs )
        derivs[p] = (xlylzl1 - calc) / step
        dargs[p] = psave
    return calc, derivs


oripnames = [ 'wedge', 'chi' ] 
oripsteps = {'wedge' : lambda wedge : 0.1,
              'chi' : lambda chi : 0.1,
              't_x' : lambda t : 0.1, 
              't_y' : lambda t : 0.1, 
              't_z' : lambda t : 0.1, 
             }

def grainorigins( colfile, agrain, pars ):
    omega = colfile.omega
    dargs = { 'wedge' : pars.get('wedge'),
                'chi' : pars.get('chi') }
    dargs['t_x'], dargs['t_y'], dargs['t_z'] = \
                agrain.translation
    origins = ImageD11.transform.compute_grain_origins( 
                omega, **dargs )
    return origins

def derivativeorigins(colfile, agrain, pars, refinables):
    origins0 = grainorigins( colfile, agrain, pars )
    omega = colfile.omega
    t_x, t_y, t_z =  agrain.translation
    wedge = pars.get('wedge')
    chi = pars.get('chi')
    derivs = {}
    dargs = { 'wedge' : pars.get('wedge'),
                'chi' : pars.get('chi') }
    dargs['t_x'], dargs['t_y'], dargs['t_z'] = \
                agrain.translation
    for p in refinables:
        step = oripsteps[p](dargs[p])
        psave = dargs[p]
        dargs[p] += step
        origins1 = ImageD11.transform.compute_grain_origins( omega, 
                 **dargs )
        derivs[p] = (origins1 - origins0) / step
        dargs[p] = psave
    return origins0, derivs

def labvec( colfile, agrain, pars):
    v_xlylzl = computeXLYLZL( colfile, pars )
    v_xgygzg = grainorigins( colfile, agrain, pars )
    v_labvec = v_xlylzl - v_xgygzg
    return v_labvec

def derivativelabvec( colfile, agrain, pars, refinables):
    dxlp = filter( lambda p: p in detpsteps, refinables )
    xlylzl, Dxlylzl = derivativeXLYLZL( colfile, pars, dxlp )
    dxgp = filter( lambda p: p in oripsteps, refinables )
    xgygzg, Dxgygzg = derivativeorigins( colfile, agrain, pars, dxgp)
    labvec = xlylzl - xgygzg
    Dlabvec = {}
    for p in dxlp:
        Dlabvec[p] = Dxlylzl[p]
    for p in dxgp:
        if Dlabvec.has_key(p): # not normally, but for forms sake
            Dlabvec[p] -= Dxgygzg[p]
        else:
            Dlabvec[p] = -Dxgygzg[p]
    return labvec, Dlabvec

def DgDk( colfile, pars):
    omega = colfile.omega
    # this would applies a rotation matrix which depends on omega, wedge, chi
    wedge = pars.get('wedge')
    chi = pars.get('chi')
    dummy = np.zeros( (3, len(omega) ) )
    dummy[0,:]=1
    dgdk0 = ImageD11.transform.compute_g_from_k(dummy,
            omega, wedge=wedge, chi=chi )
    dummy[0,:]=0
    dummy[1,:]=1
    dgdk1 = ImageD11.transform.compute_g_from_k(dummy,
            omega, wedge=wedge, chi=chi )
    dummy[1,:]=0
    dummy[2,:]=1
    dgdk2 = ImageD11.transform.compute_g_from_k(dummy,
            omega, wedge=wedge, chi=chi )
    return dgdk0,dgdk1,dgdk2


def g_from_k( kvecs, colfile, pars ):
    d0,d1,d2 = DgDk( colfile, pars ) 
    g = kvecs[0,:]*d0 + kvecs[1,:]*d1 + kvecs[2,:]*d2
    return g

def derivative_g_from_k( kvecs, Dkvecs, colfile, pars, refinables ):
    omega = colfile.omega
    # this would applies a rotation matrix which depends on omega, wedge, chi
    wedge = pars.get('wedge')
    chi = pars.get('chi')
    dummy = np.zeros( (3, len(omega) ) )
    dummy[0,:]=1
    dgdk0 = ImageD11.transform.compute_g_from_k(dummy,
            omega, wedge=wedge, chi=chi )
    dummy[0,:]=0
    dummy[1,:]=1
    dgdk1 = ImageD11.transform.compute_g_from_k(dummy,
            omega, wedge=wedge, chi=chi )
    dummy[1,:]=0
    dummy[2,:]=1
    dgdk2 = ImageD11.transform.compute_g_from_k(dummy,
            omega, wedge=wedge, chi=chi )
    g = kvecs[0,:]*dgdk0 + kvecs[1,:]*dgdk1 + kvecs[2,:]*dgdk2
    Dg = {}
    if 'wedge' in refinables:
        dummy[2,:]=0
        dummy[0,:]=1
        sw= oripsteps['wedge'](wedge)
        mat00 = ImageD11.transform.compute_g_from_k(dummy,
            omega, wedge=wedge+sw, chi=chi )
        dummy[0,:]=0
        dummy[1,:]=1
        mat10 = ImageD11.transform.compute_g_from_k(dummy,
            omega, wedge=wedge+sw, chi=chi )
        dummy[1,:]=0
        dummy[2,:]=1
        mat20 = ImageD11.transform.compute_g_from_k(dummy,
            omega, wedge=wedge+sw, chi=chi )
        dmat0dw = (mat00 - dgdk0)/sw
        dmat1dw = (mat10 - dgdk1)/sw
        dmat2dw = (mat20 - dgdk2)/sw
    if 'chi' in refinables:
        dummy[2,:]=0
        dummy[0,:]=1
        sc = oripsteps['chi'](chi)
        mat01 = ImageD11.transform.compute_g_from_k(dummy,
            omega, wedge=wedge, chi=chi+sc )
        dummy[0,:]=0
        dummy[1,:]=1
        mat11 = ImageD11.transform.compute_g_from_k(dummy,
            omega, wedge=wedge, chi=chi+sc )
        dummy[1,:]=0
        dummy[2,:]=1
        mat21 = ImageD11.transform.compute_g_from_k(dummy,
            omega, wedge=wedge, chi=chi+sc )
        dmat0dc = (mat01 - dgdk0)/sc
        dmat1dc = (mat11 - dgdk1)/sc
        dmat2dc = (mat21 - dgdk2)/sc
    for p in refinables:
        dk = Dkvecs[p]
        # g = sum_i dgdki.ki
        # dgdp = sum_i dgdki/dp.ki + dgdki.dkidp
        # print p, k.shape
        Dg[p] = dk[0,:]*dgdk0 + dk[1,:]*dgdk1 + dk[2,:]*dgdk2
        if p == 'wedge':
            # add on dmat contribution
            Dg[p] += kvecs[0,:]*dmat0dw + kvecs[1,:]*dmat1dw + \
                     kvecs[2,:]*dmat2dw
        if p == 'chi':
            # add on dmat contribution
            Dg[p] += kvecs[0,:]*dmat0dc + kvecs[1,:]*dmat1dc + \
                     kvecs[2,:]*dmat2dc
    return g, Dg


def rotv3n( mat, Dmat, v3n, Dv3n ):
    """
    Rotate vector by mat
    """
    rv = np.dot( mat, v3n )
#   print rv.shape, mat.shape, v3n.shape
#   assert (rv[:,0] == np.dot(mat, v3n[:,0])).all()
    Drv = {}
    for p in Dv3n:
        Drv[p] = np.dot( mat, Dv3n[p] )
    for p in Dmat:
        assert not Drv.has_key(p)
        Drv[p] = np.dot( Dmat[p], v3n )
    return rv, Drv

def modv3n( vec3n ):
    v2 = vec3n*vec3n
    assert v2.shape[0] == 3
    modv = np.sqrt(v2.sum(axis=0)) 
    return modv


def Derivativemodv3n( vec3n, Dvec3n ):
    v2 = vec3n*vec3n
    assert v2.shape[0] == 3
    modv = np.sqrt(v2.sum(axis=0)) 
    # Now for derivatives.
    # We have vec3n and d(vec3n)/dp
    # We want modv and d(modv)/dp
    # d(modv) = sqrt( x*x+y*y+z*z )
    # d(modv)/dp = (dmodv/dx).(dx/dp) + (dmodv/dy).(dy/dp) + (dmodv/dz).(dz/dp)
    dmodvdx = vec3n / modv
    Dmodv3n={}
    for p in Dvec3n: # each parameter in turn
        Dmodv3n[p] = (dmodvdx*Dvec3n[p]).sum(axis=0)
    return modv, Dmodv3n


def quotient_v3n_v1n( vlabvec, vmodv3n ):
    quotient = vlabvec / vmodv3n 
    return quotient

def Derivative_quotient_v3n_v1n( labvec, Dlabvec, vmodv3n, Dmodv3n):
    quotient = labvec / vmodv3n 
    Dquotient = {}
    for p in Dlabvec.keys():
        Dquotient[p] = Dlabvec[p] / vmodv3n
    for p in Dmodv3n.keys():
        if Dquotient.has_key(p): 
            Dquotient[p] -= Dmodv3n[p] * labvec / vmodv3n / vmodv3n
        else:
            Dquotient[p] = -Dmodv3n[p] * labvec / vmodv3n / vmodv3n
    return quotient, Dquotient


def Derivativegobs( cf, gr, pars, refinables):
    # Diffracted ray vector, depends on pars
    # labvec is an orthogonal basis (real laboratory coordinates)
    labvec, Dlabvec = derivativelabvec( cf, gr, pars,
                           refinables )
    # Length of this vector
    valmodv3n, Dmodv3n = Derivativemodv3n( labvec, Dlabvec )
    # direction is along here: (eg, direction cosines, normalised vector)
    direction, Ddirection = Derivative_quotient_v3n_v1n( 
                        labvec, Dlabvec, valmodv3n, Dmodv3n)
    wavelength = pars.get('wavelength')
    k = direction / wavelength
    k[0,:] = k[0,:] - 1.0/wavelength # incident beam along x
    Dk = {}
    for p in Ddirection:
        Dk[p] = Ddirection[p]/wavelength
    # Now we would like to fit this in terms of the grain ubi matrix and hkl
    # indices
    # g = transform.compute_g_from_k( k , colfile.omega, pars.get("wedge"), pars.get("chi") )
    # Find hkl indices and which peaks are used
    #hklr = np.dot( grains[0].ubi, g )
    #hkli = np.floor(hklr+0.5)
    #drlv = (hklr - hkli)
    #drlv2 = np.sqrt((drlv*drlv).sum(axis=0))
    # Rotate k to g, still an orthogonal reciprocal angstrom based metric
    gobs, Dg = derivative_g_from_k( k , Dk, cf, pars,  refinables )
    return gobs, Dg

def gobs( cf, gr, pars):
    # labvec is an orthogonal basis (real laboratory coordinates)
    v_labvec = labvec( cf, gr, pars )
    valmodv3n = modv3n( v_labvec )
    # direction is along here: (eg, direction cosines, normalised vector)
    direction = quotient_v3n_v1n( v_labvec , valmodv3n )
    wavelength = pars.get('wavelength')
    k = direction / wavelength
    k[0,:] = k[0,:] - 1.0/wavelength # incident beam along x
    gobs = g_from_k( k , cf, pars )
    return gobs

def derivativegcalc(gr, gobs):
    # find hkl indices of spots
    # FIXME : these ought to be fixed from the colfile
    # gradients are in Dg
    hkli = np.floor( np.dot( gr.ubi, gobs ) + 0.5 )
    ub = np.linalg.inv( gr.ubi )
    gr.ub = ub
    gcalc = np.dot( ub, hkli )
    # Derivatives of gcalc w.r.t UB
    Dgcalc = {}
    for i in range(3):
        for j in range(3):
            t = np.zeros((3,3))
            t[i,j] = 1
            name = 'ub%d%d'%(i,j)
            Dgcalc[name] = t
    gcalc, dGcalc = rotv3n( ub, Dgcalc, hkli, {} )
    return gcalc, dGcalc
    


class tmpcol:
    def __init__(self, cf):
        self.sc = cf.sc.copy()
        self.fc = cf.fc.copy()
        self.omega = cf.omega.copy()

def dGdobs( cf, gr, pars, step=0.01) :
    tmp = tmpcol( cf )

    g0 = gobs( tmp, gr, pars )
    np.add(tmp.sc, step, tmp.sc )
    g1 =  gobs( tmp, gr, pars )
    dgdsc = ( g1 - g0 ) / step      # 3,N
    np.subtract(tmp.sc,step,tmp.sc)

    np.add(tmp.fc, step, tmp.fc )
    g1 =  gobs( tmp, gr, pars )
    dgdfc = ( g1 - g0 ) / step
    np.subtract(tmp.fc,step,tmp.fc) # 3,N

    np.add(tmp.omega, step, tmp.omega )
    g1 =  gobs( tmp, gr, pars )
    dgdomega = ( g1 - g0 ) / step   # 3,N
    
    return dgdsc, dgdfc, dgdomega


def rotate_errors( cf , gr, pars):

    # single data point, uncorrelated has a matrix of:
    #  1/ss   0     0
    #   0    1/ff   0
    #   0     0    1/oo
    #
    # correlated that would be:
    #  1/ss  1/sf  1/so
    #  1/fs  1/ff  1/fo
    #  1/os  1/of  1/oo
    weights = np.zeros( (3,3,cf.nrows ), np.float)
    weights[0,0] = 1./cf.sigs/cf.sigs
    weights[1,1] = 1./cf.sigf/cf.sigf
    weights[2,2] = 1./cf.sigo/cf.sigo
    weights[1,0] = weights[0,1]= cf.covsf/cf.sigs/cf.sigf
    weights[2,0] = weights[0,2]= cf.covso/cf.sigs/cf.sigo
    weights[2,1] = weights[1,2]= cf.covfo/cf.sigo/cf.sigf
    #
    #
    return weights



import wripaca
from ImageD11 import gv_general

def omegacalc_ub(self, ubi, h, pars, romega):
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
    return romegacalc







