#!/usr/bin/env python

from pedees.delaunay2d import Delaunay2d
from pedees.elliptic2d import Elliptic2d
from pedees.cg import Cg
from pedees.plot import Plot
import numpy

nx, ny = 3, 3
dx, dy = 1.0/float(nx), 1.0/float(ny)
pts = [ ]
for i in range(nx + 1):
  for j in range(ny + 1):
    pts.append( numpy.array([i*dx, j*dy]) )

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

print 'stiffness matrix before BCs: ', elliptic.getStiffnessMatrix()


large = 1.e6
bcs = {}
bedges = delaunay.getBoundaryEdges()
for be in bedges:
  i, j = be
  pi = delaunay.getPoints()[i]
  pj = delaunay.getPoints()[j]
  if pi[1] > 0.999 and pj[1] > 0.999:
    # Dirichlet 1
    bcs[be] = (large, large)
  elif pi[1] < 0.01 and pj[1] < 0.01:
    # Dirichlet 0
    bcs[be] = (large, 0.0)
  # zero Neumann on all other edges (no need to specify)
elliptic.applyBoundaryConditions(bcs)

mat, b = elliptic.getStiffnessMatrix(), elliptic.getSourceVector()
slvr = Cg(mat, b)
n = len(b)
p = numpy.array([mat[i, i] for i in range(n)])
zeros = numpy.zeros( b.shape, numpy.float64 )
x = slvr.solve(precond=p, x0=zeros, numIters=100, tol=1.e-10, verbose=True)
print 'solution: ', x
print 'points: ', delaunay.getPoints()
print 'BCs: ', bcs
print 'solution error: ', slvr.getSolutionError(x)

pl = Plot(delaunay, width=300, height=300)
pl.show(x)