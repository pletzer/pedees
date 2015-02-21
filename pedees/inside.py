#!/usr/bin/env python

import numpy

class Inside:

  def __init__(self, points, faces):

    self.eps = 1.23456789e-14

    self.faces = faces
    self.points = points

    # must have at least one point
    self.ndims = len(points[0])

    # matrix and right and side vector
    self.mat = numpy.zeros( (self.ndims, self.ndims), numpy.float64 )
    self.b = numpy.zeros( (self.ndims,), numpy.float64 )

    # the optimal direction for shooting the ray
    self.direction = float('inf') * numpy.ones((self.ndims,), numpy.float64)

    # find the box corners
    self.xmins = +float('inf') * numpy.ones(self.ndims, numpy.float64)
    self.xmaxs = -float('inf') * numpy.ones(self.ndims, numpy.float64)
    for i in range(self.ndims):
      xmin = min([p[i] for p in points])
      xmax = max([p[i] for p in points])
      self.xmins[i] = min(self.xmins[i], xmin)
      self.xmaxs[i] = max(self.xmaxs[i], xmax)

  def isInside(self, point):

    # quick check, point must be inside box
    outsideBox = reduce(lambda x,y: x or y, 
      [(point[i] < self.xmins[i] - self.eps) or (point[i] > self.xmaxs[i] + self.eps) 
      for i in range(self.ndims)])

    if outsideBox:
      return False

    # any direction will do but things will run faster if the direction
    # points to the nearest domain box face (fewer intersections to 
    #  compute)
    self.computeOptimalDirection(point) # self.xmaxs - point

    # set the first column in our matrix (independent of the face)
    self.mat[:, 0] = self.direction

    # compute the number of intersections between the ray and the faces
    numberOfIntersections = 0
    for face in self.faces:

      # compute the overlap betwen the ray box and the face box
      if not self.areBoxesOverlapping(point, face):
        # intersection is impossible, don't bother...
        continue

      lmbda, xis = self.computeIntersection(point, face)
      if lmbda > 0.0 + self.eps:
        # the direction is ok
        sums = numpy.cumsum(xis)
        isInside = reduce(lambda x, y: x and y, 
          [xis[i] > 0.0 + self.eps and xis[i] < 1.0 - sums[i] - self.eps \
          for i in range(len(xis))])
        if isInside:
          numberOfIntersections += 1

    return not (numberOfIntersections % 2)

  def areBoxesOverlapping(self, point, face):
    for i in range(self.ndims):

      xminFace = min( [self.points[j][i] for j in face] )
      xmaxFace = max( [self.points[j][i] for j in face] )

      p, d = point[i], self.direction[i]
      xminRay = min(p, p + d)
      xmaxRay = max(p, p + d)

      if (xmaxRay < xminFace + self.eps) or \
         (xmaxFace < xminRay + self.eps):
         # no overlap
        return False
    return True

  def computeNormalToPlane(self, axis, pm):
    normal = numpy.zeros( (self.ndims,), numpy.float64 )
    normal[axis] = pm
    loCorner = self.xmins
    if pm > 0:
      loCorner[axis] = self.xmaxs[axis]
    return normal, loCorner

  def computeOptimalDirection(self, point):
    """
    Compute the direction of the ray to to the nearest 
    domain box face. This will update self.direction
    @param point starting point of the ray
    """
    # iterate over the faces of the box then find the minimum distance between
    # the point and the face. 
    minDistance = float('inf')
    for pm in (-1, 1):
      for axis in range(self.ndims):
        normal, loCorner = self.computeNormalToPlane(axis, pm)
        distance = numpy.dot(normal, loCorner - point)
        if abs(distance) < minDistance:
          self.direction = normal * distance
          minDistance = abs(distance)

  def computeIntersection(self, point, face):
    self.b = self.points[face[0]] - point
    for i in range(self.ndims - 1):
      self.mat[:, 1 + i] = self.points[face[0]] - self.points[face[1 + i]]
    solution = numpy.linalg.solve(self.mat, self.b)
    return solution[0], solution[1:]


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


