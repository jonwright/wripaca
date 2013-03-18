

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

def anglediff(a1, a2):
    return np.arctan2(np.sin(a1-a2), np.cos(a1-a2))



for wedge, chi in [(0.,0.),(-10.,0.), (11,-15), (0,6)]:
    print "Wedge",wedge,"chi",chi
    print "gve[0] gve[1] gve[2]  om1py  om2py  om1af  om2af  omtrue"
    gve = transform.compute_g_vectors( tth, eta, omega, wvln, wedge, chi )
    gve = gve.T.copy()
    post = gv_general.wedgechi( wedge, chi ) 
    posti = np.linalg.inv( post )

    # don't try to match this as it is wrong!!
    #    tth, eta, omega = transform.uncompute_g_vectors( hkl.T, wvln, 
    #        wedge=wedge,
    #        chi = chi)

    ppre = gv_general.axis_from_matrix(np.reshape(pre,(3,3)))
    ppost = gv_general.axis_from_matrix(np.reshape(post,(3,3)))


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

