

import sys
sys.path.append("build/lib.win-amd64-2.7")

import wripaca
print wripaca.__file__

import numpy as np
from ImageD11 import transform, gv_general

print "det(eye(3))",wripaca.determinant3( np.eye(3).ravel() )

pre = np.eye(3).ravel()
post = np.eye(3).ravel()
axis =  np.array( (0.,0.,-1.))

tth = np.array( [60, 10, 15, 20,25, 1.9555],np.float)
wvln = 0.5
eta = np.array( [90, 10,120,-20,340, -73 ],np.float)
omega = np.array([60,90,180, 60,97,  131],np.float)

def testaxisang():
    v = np.array([11.,12.,13.],np.float)
    cv = v.copy()
    axis = np.array([np.sqrt(1./3),-np.sqrt(1./3),np.sqrt(1./3)],np.float)

    for i in range(len(omega)):
        a = gv_general.rotation_axis(axis, angle=omega[i])
        rv = a.rotate_vectors( v )
        wripaca.rotate_vector_axis_angle(
                axis,
                omega[i]*np.pi/180,
                v,
                cv)
        co=np.cos(omega[i]*np.pi/180)
        so=np.sin(omega[i]*np.pi/180)
        M = np.array([[co,so,0],[-so,co,0],[0,0,1]])
        sv= np.dot(M,v)
        if 0:
            print "Omega =",omega[i]
            print "py    ",rv
            print "simple",sv
            print "affine",cv
        assert ((cv-rv)**2).sum()<1e-10,"badness in axis angle"
    print "AxisAngleOK"

testaxisang()



def anglediff(a1, a2):
    return np.arctan2(np.sin(a1-a2), np.cos(a1-a2))



for wedge, chi in [(0.,0.),(-10.,0.), (11,-15), (0,6)]:
    print "Wedge",wedge,"chi",chi
    print "gve[0] gve[1] gve[2]  om1py  om2py  om1af  om2af  omtrue"
    gve = transform.compute_g_vectors( tth, eta, omega, wvln, wedge, chi )
    gve = gve.T.copy()
    post = gv_general.chiwedge( wedge=wedge, chi=chi ) 
    posti = np.dot(gv_general.wedgemat(wedge), gv_general.chimat(chi)).T

    # don't try to match this as it is wrong!!
    #    tth, eta, omega = transform.uncompute_g_vectors( hkl.T, wvln, 
    #        wedge=wedge,
    #        chi = chi)

    ppre = np.reshape(pre,(3,3))
    ppost = np.reshape(post,(3,3))


    omg1, omg2, valid = gv_general.g_to_k( gve.T, wvln, axis=axis, pre=ppre, post=ppost )
    #print pre
    #print posti.ravel()
    #print axis
    #print wvln
    for i in range(len(gve)):
        print "%8.4f%8.4f%8.4f"%tuple(gve[i]),
        print "%9.3f %9.3f"%(omg1[i], omg2[i]), 
        #print gve.T[i],
        om1, om2 = wripaca.omegacalc( gve[i], pre, posti.ravel(), axis, wvln )
        print "%9.3f %9.3f"%(om1*180/np.pi, om2*180/np.pi), omega[i]

        error =min( abs(anglediff( omega[i]*np.pi/180, om1 )), \
                    abs(anglediff( omega[i]*np.pi/180, om2 ) ))
        assert error < 1e-15
#        print error
    print 



def testhklcalc_many(wedge, chi):
    # hklcalc_many
    # apply wedge, chi to axis if necessary

    axis = np.array( [0,0,-1], np.float )
    UBI = np.eye(3).ravel()
    hkl = np.zeros( (len(omega), 3), np.float)
    XL = np.ones( (len(omega), 3), np.float)*1e3
    T = np.array([100,200,300],np.float)
    wvln = 0.254
    pre = np.eye(3).ravel()
    post = gv_general.chiwedge( wedge=wedge, chi=chi ).ravel()
    wripaca.hklcalc_many( XL, axis,  omega*np.pi/180, UBI, T, pre, post, hkl, wvln )
    print "axis is",axis
    t, e = transform.compute_tth_eta_from_xyz( XL.T,
        t_x = T[0],
        t_y = T[1],
        t_z = T[2],
        omega = omega,
        wedge = wedge,
        chi = chi
        )
    #print t,e,omega
    ori = transform.compute_grain_origins( omega, wedge=wedge, chi=chi, 
        t_x  = T[0], t_y = T[1], t_z = T[2] )
    print "last origin",ori.T[-1],ori.T[0]
    print "last xl",XL[-1],XL[0]
    kve = transform.compute_k_vectors( t, e, wvln )
    print "last k",kve.T[-1],kve.T[0]
    gve = transform.compute_g_vectors( t, e, omega, wvln, wedge=wedge, chi=chi )
    print "last g",gve.T[-1],gve.T[0]
    htest = np.dot( UBI.reshape(3,3), gve )
    for i in range(len(t)):
        #        print gve[:,i]
        # print hkl[i]

        if ((gve[:,i]-hkl[i])**2).sum() > 1e-10:
            print "FAIL",
            print i,gve[:,i],hkl[i],((gve[:,i]-hkl[i])**2).sum() 
            print wedge,chi
            sys.exit()
    
    print "OK for wedge",wedge,"chi",chi


testhklcalc_many( 10.0, 0.0 )
testhklcalc_many( 0.0, 0.0 )
testhklcalc_many( 0.0, 10.0 )
testhklcalc_many( 10.0, 11.0)

