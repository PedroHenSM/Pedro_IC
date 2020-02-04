# File : setup.py 
  
from distutils.core import setup, Extension 
#name of module 
name  = "cec20NoConstraints"
  
#version of module 
version = "1.0"
  
# specify the name of the extension and source files 
# required to compile this 
# ext_modules = Extension(name='_cec20NoConstraints',sources=["_cec20NoConstraints_module.cc","C version/cec20_test_func.cpp", "main.cpp"]) 
ext_modules = Extension(name='_cec20NoConstraints',sources=["_cec20NoConstraints_module.cc", "C version/cec20_test_func.cpp", "C version/main.cpp"]) 
  
setup(name=name, 
      version=version, 
      ext_modules=[ext_modules])


"""
***COMPILING***
swig -python -c++ -o _cec20NoConstraints_module.cc cec17NoConstraints.i
python3 setupcec17NoConstraints.py build_ext --inplace

***IMPORTING***
import sys
sys.path.insert(0,'nameOfFolder')
import nameOfFile

"""