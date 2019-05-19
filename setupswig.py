# File : setup.py 
# Create the example.i before using this script  
from distutils.core import setup, Extension 
#name of module 
name  = "example"
  
#version of module 
version = "1.0"
  
# specify the name of the extension and source files 
# required to compile this 
ext_modules = Extension(name='_example',sources=["_example_module.cc","example1.cpp","exemple2.cpp","soOn.cpp"]) 
  
setup(name=name, 
      version=version, 
      ext_modules=[ext_modules]) 


"""
Note: on "swig -python" must be python2
Note: on "python setupexample.py" must be python3 (if you are using python3)
***COMPILING***
swig -python -c++ -o _example_module.cc example.i
python setupexample.py build_ext --inplace

***IMPORTING***
import sys
sys.path.insert(0,'nameOfFolder')
import nameOfFile

"""