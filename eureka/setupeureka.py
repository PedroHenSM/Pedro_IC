# File : setup.py 
  
from distutils.core import setup, Extension 
#name of module 
name  = "eureka"
  
#version of module 
version = "1.0"
  
# specify the name of the extension and source files 
# required to compile this 
ext_modules = Extension(name='_eureka',sources=["_eureka_module.cc","EurekaOptimaException.cpp","Problem.cpp","TrussBarStructureStaticProblem.cpp","TrussBarStructureStaticSimulator.cpp","F101Truss10Bar.cpp","F103Truss25Bar.cpp","F105Truss60Bar.cpp","F107Truss72Bar.cpp","F109Truss942Bar.cpp"]) 
  
setup(name=name, 
      version=version, 
      ext_modules=[ext_modules]) 


"""
***COMPILING***
swig -python -c++ -o _eureka_module.cc eureka.i
python setupeureka.py build_ext --inplace

***IMPORTING***
import sys
sys.path.insert(0,'nameOfFolder')
import nameOfFile

"""