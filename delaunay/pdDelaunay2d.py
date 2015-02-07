#!/usr/bin/env python

import numpy
import math
import copy

class Delaunay2d:

  EPS = 1.23456789e-14

  def __init__(self, points, boundaryEdges=None, holes=None):

    # data structures
    self.points = points[:] # copy
    self.triangles = [] # cells
    self.edge2Triangles = {} # edge to triangle(s) map
    self.boundaryEdges = set()
    self.holes = holes
    
    # compute center of gravity
    cg = numpy.zeros((2,), numpy.float64)
    for pt in points:
      cg += pt
    cg /= len(points)

    # sort
    def distanceSquare(pt):
      d = pt - cg
      return numpy.dot(d, d)
    self.points.sort(key = distanceSquare)

    # create first triangle, make sure we're getting a non-zero area otherwise
    # drop the points
    area = 0.0
    index = 0
    stop = False
    while not stop and index + 2 < len(points):
      area = self.getArea(index, index + 1, index + 2)
      if abs(area) < self.EPS:
        del self.points[index]
      else:
        stop = True
    if index <= len(self.points) - 3:
      tri = [index, index + 1, index + 2]
      self.makeCounterClockwise(tri)
      self.triangles.append(tri)
      # boundary edges
      e01 = (0, 1)
      self.boundaryEdges.add(e01)
      e12 = (1, 2)
      self.boundaryEdges.add(e12)
      e20 = (2, 0)
      self.boundaryEdges.add(e20)
      self.edge2Triangles[e01] = [0,]
      self.edge2Triangles[e12] = [0,]
      e02 = (0, 2)
      self.edge2Triangles[e02] = [0,]

    else:
      # all the points fall on a line
      return

    # add additional points
    for i in range(3, len(self.points)):
      self.addPoint(i)

    # remove all triangles inside holes
    # TO DO 

  def getTriangles(self):
    return self.triangles

  def getEdges(self):
    return self.edge2Triangles.keys()

  def getArea(self, ip0, ip1, ip2):
    """
    Compute the parallelipiped area
    @param ip0 index of first vertex
    @param ip1 index of second vertex
    @param ip2 index of third vertex
    """
    d1 = self.points[ip1] - self.points[ip0]
    d2 = self.points[ip2] - self.points[ip0]
    return (d1[0]*d2[1] - d1[1]*d2[0])

  def isEdgeVisible(self, ip, edge):
    """
    Return true iff the point lies to its right when the edge points down
    @param ip point index
    @param edge (2 point indices with orientation)
    @return True if visible    
    """
    area = self.getArea(ip, edge[0], edge[1])
    if area < self.EPS:
      return True
    return False

  def makeCounterClockwise(self, ips):
    area = self.getArea(ips[0], ips[1], ips[2])
    if area < -self.EPS:
      ip1, ip2 = ips[1], ips[2]
      # swap
      ips[1], ips[2] = ip2, ip1

  def flipEdges(self):
    res = False
    edgesToRemove = []
    edgesToAdd = {}
    for edge, tris in self.edge2Triangles.items():
      if len(tris) < 2:
        continue
      iTri1, iTri2 = tris
      tri1 = self.triangles[iTri1]
      tri2 = self.triangles[iTri2]
      iOpposite1 = -1
      iOpposite2 = -1
      for i in range(3):
        if not tri1[i] in edge:
          iOpposite1 = tri1[i]
        if not tri2[i] in edge:
          iOpposite2 = tri2[i]
      da1 = self.points[edge[0]] - self.points[iOpposite1]
      db1 = self.points[edge[1]] - self.points[iOpposite1]
      da2 = self.points[edge[0]] - self.points[iOpposite2]
      db2 = self.points[edge[1]] - self.points[iOpposite2]
      crossProd1 = self.getArea(iOpposite1, edge[0], edge[1])
      crossProd2 = self.getArea(iOpposite2, edge[1], edge[0])
      dotProd1 = numpy.dot(da1, db1)
      dotProd2 = numpy.dot(da2, db2)
      angle1 = abs(math.atan2(crossProd1, dotProd1))
      angle2 = abs(math.atan2(crossProd2, dotProd2))
      if angle1 + angle2 > math.pi*(1.0 + self.EPS):
        # flip the triangles
        newTri1 = [iOpposite1, edge[0], iOpposite2]
        newTri2 = [iOpposite1, iOpposite2, edge[1]]
        self.triangles[iTri1] = newTri1
        self.triangles[iTri2] = newTri2
        edgesToRemove.append(edge)
        e = tuple([iOpposite1, iOpposite2].sort())
        edgesToAdd[e] = [newTri1, newTri2]
        res = True
      # remove edges
      for e in edgesToRemove:
        del self.edge2Triangles[e]
      # add edges
      for e, tris in edgesToAdd.items():
        self.edge2Triangles[e] = tris
    return res

  def addPoint(self, ip):

    for edge in copy.copy(self.boundaryEdges):

      if self.isEdgeVisible(ip, edge):

        # create new triangle
        newTri = [edge[0], edge[1], ip]
        newTri.sort()
        self.makeCounterClockwise(newTri)
        self.triangles.append(newTri)

        # update the edge to triangle map
        e = list(edge[:])
        e.sort()
        iTri = len(self.triangles) - 1 
        self.edge2Triangles[tuple(e)].append(iTri)

        # add the two boundary edges
        e1 = [ip, edge[0]]
        e1.sort()
        e1 = tuple(e1)
        e2 = [edge[1], ip]
        e2.sort()
        e2 = tuple(e2)
        self.edge2Triangles[e1] = [iTri,]
        self.edge2Triangles[e2] = [iTri,]

        # update the boundary edges
        self.boundaryEdges.remove(edge)
        self.boundaryEdges.add(e1)
        self.boundaryEdges.add(e2)

    # recursively flip edges
    flipped = True
    while flipped:
      flipped = self.flipEdges()

#############################################################################

def testOneTriangle():
  xyPoints = [numpy.array([0., 0.]), numpy.array([1., 0.]), numpy.array([0., 1.])]
  delaunay = Delaunay2d(xyPoints)
  print 'triangles: ', delaunay.getTriangles()
  print 'edges: ', delaunay.getEdges()

def testOneTriangle2():
  # points go clockwise
  xyPoints = [numpy.array([0., 0.]), numpy.array([0., 1.]), numpy.array([1., 0.])]
  delaunay = Delaunay2d(xyPoints)
  print 'triangles: ', delaunay.getTriangles()
  print 'edges: ', delaunay.getEdges()

def testTwoTriangles():
  xyPoints = [numpy.array([0., 0.]), numpy.array([1., 0.]), numpy.array([0., 1.]), numpy.array([1., 1.])]
  delaunay = Delaunay2d(xyPoints)
  print 'triangles: ', delaunay.getTriangles()
  print 'edges: ', delaunay.getEdges()


if __name__ == '__main__': 
  #testOneTriangle()
  #testOneTriangle2()
  testTwoTriangles()



