

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from numpy.distutils.misc_util import get_numpy_include_dirs

nid = get_numpy_include_dirs() + ["../include"]
print nid
e = Extension("wripaca", 
            ["wripaca.pyx", "cwripaca.pxd",
            "../src/affine.c"],# , "../src/wripaca.c"],
        include_dirs = nid ,
        )


setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [e]
)



