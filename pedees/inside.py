#!/usr/bin/env python

import numpy

class Inside:

  def __init__(self, points, faces):

    self.eps = 1.23456789e-14

    self.faces = faces
    self.points = points

    # must have at least one point
    self.ndims = len(points[0])

    self.mat = numpy.zeros( (self.ndims, self.ndims), numpy.float64 )
    self.b = numpy.zeros( (self.ndims,), numpy.float64 )

    # find the box corners
    self.xmins = +float('inf') * numpy.ones(self.ndims, numpy.float64)
    self.xmaxs = -float('inf') * numpy.ones(self.ndims, numpy.float64)
    for i in range(self.ndims):
      xmin = min([p[i] for p in points])
      xmax = max([p[i] for p in points])
      self.xmins[i] = min(self.xmins[i], xmin)
      self.xmaxs[i] = max(self.xmaxs[i], xmax)

  def isPointRoughlyInFaceBox(self, point, face):
    for i in range(self.ndims):
      xmin = min( [self.points[j][i] for j in face] )
      if point[i] < xmin - self.eps:
        return False
      xmax = max( [self.points[j][i] for j in face] )
      if point[i] > xmax + self.eps:
        return False
    return True

  def computeIntersection(self, point, face, direction):
    self.mat[:, 0] = direction
    self.b = self.points[face[0]] - point
    for i in range(self.ndims - 1):
      self.mat[:, 1 + i] = self.points[face[0]] - self.points[face[1 + i]]
    solution = numpy.linalg.solve(self.mat, self.b)
    return solution[0], solution[1:]

  def isInside(self, point):

    # quick check first, point must be inside box
    outsideBox = reduce(lambda x,y: x or y, 
      [(point[i] < self.xmins[i] - self.eps) or (point[i] > self.xmaxs[i] + self.eps) 
      for i in range(self.ndims)])

    if outsideBox:
      return False

    # direction is towards the high end box corner
    direction = self.xmaxs - point
    numberOfIntersections = 0
    for face in self.faces:

      if not self.isPointRoughlyInFaceBox(point, face):
        continue

      lmbda, xis = self.computeIntersection(point, face, direction)
      if lmbda > 0.0 + self.eps:
        # the direction is right
        sums = numpy.cumsum(xis)
        isInside = reduce(lambda x, y: x and y, 
          [xis[i] > 0.0 + self.eps and xis[i] < 1.0 - sums[i] - self.eps \
          for i in range(len(xis))])
        if isInside:
          numberOfIntersections += 1

    return not (numberOfIntersections % 2)

########################################################################################

def test2d_1():
  points = [numpy.array([0., 0.]),
            numpy.array([1., 0.]),
            numpy.array([0., 1.])]
  faces = [(0, 1), (1, 2), (2, 0)]
  inside = Inside(points, faces)

  # really outside
  assert(inside.isInside(numpy.array([-0.1, -0.2])) == False)

  # inside
  assert(inside.isInside(numpy.array([+0.1, +0.2])) == True)

  # outside
  assert(inside.isInside(numpy.array([-0.1, +0.2])) == False)
  assert(inside.isInside(numpy.array([+0.1, -0.2])) == False)

  # on the face
  assert(inside.isInside(numpy.array([-0.000001, +0.2])) == False)

if __name__ == '__main__':
  test2d_1()

