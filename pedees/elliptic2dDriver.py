#/usr/bin/env python

from delaunay2d import Delaunay2d
from inside import Inside
from elliptic2d import Elliptic2d
from cg import Cg
import math
import numpy

class Elliptic2DDriver:

  def __init__(self, fFunc, gFunc, sFunc, xFunc, yFunc):
    self.fFunc = fFunc
    self.gFunc = gFunc
    self.sFunc = sFunc
    self.xFunc = xFunc
    self.yFunc = yFunc

    self.delny = None
    self.ellipt = Elliptic2d(self.fFunc, self.gFunc, self.sFunc)
    self.slvr = None
    self.solution = None

  def triangulate(self, numCells=100):

    # number of boundary points
    self.nt = 2*int(math.sqrt(numCells))
    self.dt = 1.0/float(self.nt)

    # compute the boundary points
    points = []
    xmins =  float('inf') * numpy.ones( (2,), numpy.float64 )
    xmaxs = -float('inf') * numpy.ones( (2,), numpy.float64 )
    for i in range(self.nt):
      t = i * self.dt
      points.append( numpy.array([self.xFunc(t), self.yFunc(t)]) )

    # rough estimate of the max triangle area
    maxArea = (xmaxs[0] - xmins[0])*(xmaxs[1] - xmins[1])/float(numCells)
    self.delny = Delaunay2d(points, maxArea=maxArea)
    self.delny.triangulate()
    self.delny.makeDelaunay()
    self.delny.refine()
    self.delny.makeDelaunay() # this may not be necessary

    self.ellipt.assemble(self.delny)

  def ffunc(self, xy):
    x, y = xy
    return eval(self.fFunc)

  def gfunc(self, xy):
    x, y = xy
    return eval(self.gFunc)

  def sfunc(self, xy):
    x, y = xy
    return eval(self.sFunc)

  def applyBoundaryConditions(self, bFunc, cFunc):

    bcs = {}
    for i in range(self.nt):
      t = (i + 0.5) * self.dt
      b, c = bFunc(t), cFunc(t)
      bcs[i, (i+1)%self.nt] = (b, c)

    self.ellipt.applyBoundaryConditions(bcs)

    self.slvr = Cg(self.ellipt.getStiffnessMatrix(), \
                              self.ellipt.getSourceVector())


  def solve(self):
        
    # initial guess
    self.n = len(self.delny.points)

    mat = self.ellipt.getStiffnessMatrix()

    # use diagonal matrix as preconditioner
    precond = numpy.array( [mat[i, i] for i in range(self.n)] )

    # initial guess
    x0 = numpy.zeros( (self.n,), numpy.float64 )

    # max number of iterations
    numIters = self.n

    self.solution = self.slvr.solve(precond = precond, x0 = x0, 
                           numIters=numIters, tol=1.e-10, verbose=True)


  def show(self):
		from plot import Plot
		pl = Plot(self.delny, width=500, height=500)
		pl.show(self.solution)
    
      
################################################################################
def main():

  def fFunc(p):
    return 1.0

  def gFunc(p):
    return 0.0

  def sFunc(p):
    return 0.0

  def xFunc(t):
    return math.cos(2 * math.pi * t)

  def yFunc(t):
    return math.sin(2 * math.pi * t)

  def bFunc(t):
    return 1.0

  def cFunc(t):
    return math.sin(2 * math.pi * t)

  e = Elliptic2DDriver(fFunc, gFunc, sFunc, xFunc, yFunc)
  e.triangulate(numCells = 100)

  e.applyBoundaryConditions(bFunc, cFunc)

  e.solve()

  e.show()

if __name__ == '__main__':
  main()
