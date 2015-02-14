#!/usr/bin/env python

from pedees.delaunay2d import Delaunay2d
from pedees.elliptic2d import Elliptic2d
from pedees.cg import Cg
from pedees.plot import Plot
import numpy

nx, ny = 2, 1
dx, dy = 1.0/float(nx), 1.0/float(ny)
pts = [ numpy.array([i*dx, 0.0]) for i in range(nx) ] + \
      [ numpy.array([1.0, j*dy]) for j in range(ny) ] + \
      [ numpy.array([1.0 - i*dx, 1.0]) for i in range(nx) ] + \
      [ numpy.array([0.0, 1.0 - j*dy]) for j in range(ny) ]
delaunay = Delaunay2d(pts)
print delaunay.getPoints()
print delaunay.getTriangles()

def fFunc(xy):
	return 1.0

def gFunc(xy):
	return 0.0

def sFunc(xy):
	return 0.0

elliptic = Elliptic2d(fFunc, gFunc, sFunc)
elliptic.assemble(delaunay)

large = 1.e6
bcs = {}
bedges = delaunay.getBoundaryEdges()
for be in bedges:
  i, j = be
  pi = delaunay.getPoints()[i]
  pj = delaunay.getPoints()[j]
  if pi[0] > 0.999 and pj[0] > 0.999:
    # Dirichlet 1
    bcs[be] = (large, large)
  elif pi[0] < 0.01 and pj[0] < 0.01:
    # Dirichlet 0
    bcs[be] = (large, 0.0)
  # zero Neumann on all other edges (no need to specify)
elliptic.applyBoundaryConditions(bcs)

mat, b = elliptic.getStiffnessMatrix(), elliptic.getSourceVector()
slvr = Cg(mat, b)
n = len(b)
p = numpy.array([mat[i, i] for i in range(n)])
zeros = numpy.zeros( b.shape, numpy.float64 )
x = slvr.solve(precond=p, x0=zeros, numIters=10, tol=1.e-10, verbose=True)
print x

pl = Plot(delaunay, width=500, height=300)
pl.show(x)