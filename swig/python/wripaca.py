# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.
# This file is compatible with both classic and new-style classes.

from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_wripaca', [dirname(__file__)])
        except ImportError:
            import _wripaca
            return _wripaca
        if fp is not None:
            try:
                _mod = imp.load_module('_wripaca', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _wripaca = swig_import_helper()
    del swig_import_helper
else:
    import _wripaca
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


__author__="Jon Wright <wright@esrf.fr>"
__date__="March 2013"
__version__="?"
__doc__="..."




def new_doubleArray(*args):
  return _wripaca.new_doubleArray(*args)
new_doubleArray = _wripaca.new_doubleArray

def delete_doubleArray(*args):
  return _wripaca.delete_doubleArray(*args)
delete_doubleArray = _wripaca.delete_doubleArray

def doubleArray_getitem(*args):
  return _wripaca.doubleArray_getitem(*args)
doubleArray_getitem = _wripaca.doubleArray_getitem

def doubleArray_setitem(*args):
  return _wripaca.doubleArray_setitem(*args)
doubleArray_setitem = _wripaca.doubleArray_setitem
class afmatrix(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, afmatrix, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, afmatrix, name)
    __repr__ = _swig_repr
    __swig_setmethods__["m"] = _wripaca.afmatrix_m_set
    __swig_getmethods__["m"] = _wripaca.afmatrix_m_get
    if _newclass:m = _swig_property(_wripaca.afmatrix_m_get, _wripaca.afmatrix_m_set)
    __swig_setmethods__["v"] = _wripaca.afmatrix_v_set
    __swig_getmethods__["v"] = _wripaca.afmatrix_v_get
    if _newclass:v = _swig_property(_wripaca.afmatrix_v_get, _wripaca.afmatrix_v_set)
    def __init__(self): 
        this = _wripaca.new_afmatrix()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _wripaca.delete_afmatrix
    __del__ = lambda self : None;
afmatrix_swigregister = _wripaca.afmatrix_swigregister
afmatrix_swigregister(afmatrix)


def determinant3(*args):
  return _wripaca.determinant3(*args)
determinant3 = _wripaca.determinant3


