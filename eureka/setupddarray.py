# File : setup.py 
  
from distutils.core import setup, Extension 
#name of module 
name  = "ddarray"
  
#version of module 
version = "1.0"
  
# specify the name of the extension and source files 
# required to compile this 
ext_modules = Extension(name='_ddarray',sources=["_ddarray_module.cc","ddarray.cpp"]) 
  
setup(name=name, 
      version=version, 
      ext_modules=[ext_modules]) 


"""
***COMPILING***
swig -python -c++ -o _ddarray_module.cc ddarray.i
python setupddarray.py build_ext --inplace

***IMPORTING***
import sys
sys.path.insert(0,'nameOfFolder')
import nameOfFile

"""