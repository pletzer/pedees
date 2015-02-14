#!/usr/bin/python 

from distutils.core import setup

setup(name='pedees',
      version='0.1',
      description='A simple PDE solver',
      author='Alexander Pletzer',
      author_email='alexander@gokliya.net',
      packages=['delaunay', 'fe', 'solver'],
     )