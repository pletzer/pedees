#!/usr/bin/python

from distutils.core import setup

setup(name='pedees',
      version='0.1',
      description='A simple PDE solver',
      author='Alexander Pletzer',
      author_email='alexander@gokliya.net',
      py_modules=['pedees.delaunay', 
                  'pedees.fe', 
                  'pedees.solver'],
      #package_dir={'': '../pedees'},
     )