# File : setup.py 
  
from distutils.core import setup, Extension 
#name of module 
name  = "cec17NoConstraints"
  
#version of module 
version = "1.0"
  
# specify the name of the extension and source files 
# required to compile this 
ext_modules = Extension(name='_cec17NoConstraints',sources=["_cec17NoConstraints_module.cc","cec17_test_func.cpp", "main.cpp"]) 
  
setup(name=name, 
      version=version, 
      ext_modules=[ext_modules]) 


"""
***COMPILING***
swig -python -c++ -o _cec17NoConstraints_module.cc cec17NoConstraints.i
python setupcec17NoConstraints.py build_ext --inplace

***IMPORTING***
import sys
sys.path.insert(0,'nameOfFolder')
import nameOfFile

"""