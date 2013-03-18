#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext
copt =  {    'msvc': ['/openmp', '/Ox', '/fp:fast','/favor:INTEL64','/Og']  ,
         'mingw32' : ['-fopenmp','-O3','-ffast-math','-march=native']       }
lopt =  {'mingw32' : ['-fopenmp'] }

class build_ext_subclass( build_ext ):
    def build_extensions(self):
        c = self.compiler.compiler_type
        if copt.has_key(c):
            for e in self.extensions:
                e.extra_compile_args = copt[ c ]
        if lopt.has_key(c):
            for e in self.extensions:
                e.extra_link_args = lopt[ c ]
        print self.extensions[0].extra_compile_args
        build_ext.build_extensions(self)

mod = Extension('_wripaca',
                sources=['../wripaca_wrap.c', 
                         '../../src/wripaca.c',
                         '../../src/affine.c'],
                include_dirs=['../../include']
                )

setup (name = 'wripaca',
       ext_modules = [mod],
       py_modules = ["wripaca"],
       cmdclass = {'build_ext': build_ext_subclass }
       )

