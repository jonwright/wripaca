

import sys
sys.path.append("build/lib.win-amd64-2.7")

import wripaca
print wripaca.__file__

import numpy as np

print "det(eye(3))",wripaca.determinant3( np.eye(3).ravel() )

hkl = np.array((1.,1.,1.))
UBI = np.eye(3).ravel()
pre = np.eye(3).ravel()
post = np.eye(3).ravel()
axis =  np.array( (0.,0.,1.))
wvln = 0.1

print wripaca.omegacalc( hkl, UBI, pre, post, axis, wvln )

