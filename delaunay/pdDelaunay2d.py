#!/usr/bin/env python

import numpy
import math
import copy

class Delaunay2d:

  EPS = 1.23456789e-14

  def __init__(self, points):

    # data structures
    self.points = points[:] # copy
    self.triangles = [] # cells
    self.edge2Triangles = {} # edge to triangle(s) map
    self.boundaryEdges = set()
    self.appliedBoundaryEdges = None
    self.holes = None

    
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
      e01 = (tri[0], tri[1])
      self.boundaryEdges.add(e01)
      e12 = (tri[1], tri[2])
      self.boundaryEdges.add(e12)
      e20 = (tri[2], tri[0])
      self.boundaryEdges.add(e20)
      e01 = list(e01)
      e01.sort()
      e01 = tuple(e01)
      self.edge2Triangles[e01] = [0,]
      e12 = list(e12)
      e12.sort()
      e12 = tuple(e12)
      self.edge2Triangles[e12] = [0,]
      e20 = list(e20)
      e20.sort()
      e20 = tuple(e20)
      self.edge2Triangles[e20] = [0,]

    else:
      # all the points fall on a line
      return

    print '.... after 1st tri: ', self.triangles, self.edge2Triangles, self.boundaryEdges
    # add additional points
    for i in range(3, len(self.points)):
      self.addPoint(i)
      print '.... after adding point i = ', i, self.triangles, self.edge2Triangles, self.boundaryEdges

    # remove all triangles inside holes
    # TO DO 

  def getTriangles(self):
    """
    @return triangles
    """
    return self.triangles

  def getEdges(self):
    """
    @return egdes
    """
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
    print '??????? ', ip, edge[0], edge[1], ' area = ', area
    print '??????? points: ', self.points[ip], self.points[edge[0]], self.points[edge[1]]
    if area < self.EPS:
      return True
    return False

  def makeCounterClockwise(self, ips):
    """
    Re-order nodes to ensure positive area (in-place operation)
    """
    area = self.getArea(ips[0], ips[1], ips[2])
    if area < -self.EPS:
      ip1, ip2 = ips[1], ips[2]
      # swap
      ips[1], ips[2] = ip2, ip1

  def flipEdges(self):
    """
    Flip edges to statisfy Delaunay's criterion
    """
    res = False
    return res # DEBUG

    edgesToRemove = []
    edgesToAdd = {}
    print '*** self.edge2Triangles = ', self.edge2Triangles
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
        e = [iOpposite1, iOpposite2]
        e.sort()
        e = tuple(e)
        edgesToAdd[e] = [iTri1, iTri2]
        res = True
      
      # remove edges
      print '---- edges to remove = ', edgesToRemove
      for j in range(len(edgesToRemove) - 1, -1, 1):
        e = edgesToRemove[j]
        print '---- removing ', e
        del self.edge2Triangles[e]
        del edgesToRemove[j]

      # add edges
      for e, tris in edgesToAdd.items():
        self.edge2Triangles[e] = tris
        del edgesToAdd[e]

    return res

  def addPoint(self, ip):
    """
    Add point
    @param ip point index
    """

    # collection for later updates
    boundaryEdgesToRemove = set()
    boundaryEdgesToAdd = set()

    for edge in self.boundaryEdges:

      if self.isEdgeVisible(ip, edge):

        print '>>>>>>>> point ', ip, ' sees ', edge

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
        v1 = self.edge2Triangles.get(e1, [])
        v1.append(iTri)
        v2 = self.edge2Triangles.get(e2, [])
        v2.append(iTri)
        self.edge2Triangles[e1] = v1
        self.edge2Triangles[e2] = v2

        # keep track of the boundary edges to update
        boundaryEdgesToRemove.add(edge)
        boundaryEdgesToAdd.add(e1)
        boundaryEdgesToAdd.add(e2)

    # update the boundary edges
    for bedge in boundaryEdgesToRemove:
      self.boundaryEdges.remove(bedge)
    for bedge in boundaryEdgesToAdd:
      complBEdge = list(bedge)
      complBEdge.reverse()
      complBEdge = tuple(complBEdge)
      if complBEdge not in boundaryEdgesToAdd:
        # only add boundary edge if it does not appear
        # twice in different order
        self.boundaryEdges.add(bedge)


    # recursively flip edges
    flipped = True
    while flipped:
      flipped = self.flipEdges()

  def show(self, width=500, height=300):

    import Tkinter

    xmin = min([p[0] for p in self.points])
    ymin = min([p[1] for p in self.points])
    xmax = max([p[0] for p in self.points])
    ymax = max([p[1] for p in self.points])
    w = width - 2
    h = height - 2

    master = Tkinter.Tk()
    c = Tkinter.Canvas(master, width=width, height=height)
    c.pack()
    for e in self.edge2Triangles:
      i1, i2 = e
      xp1 = 1 + int(w*(self.points[i1][0] - xmin)/(xmax - xmin))
      yp1 = 1 + int(h*(ymax - self.points[i1][1])/(ymax - ymin))
      xp2 = 1 + int(w*(self.points[i2][0] - xmin)/(xmax - xmin))
      yp2 = 1 + int(h*(ymax - self.points[i2][1])/(ymax - ymin))
      c.create_line(xp1, yp1, xp2, yp2)
    Tkinter.mainloop()

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
  #delaunay.show()

def testRandomTriangles():
  import random
  random.seed(1234)
  xyPoints = [numpy.array([random.random(), random.random()]) for i in range(5)]
  delaunay = Delaunay2d(xyPoints)
  print '*** points: ', delaunay.points
  print delaunay.edge2Triangles
  print delaunay.boundaryEdges
  delaunay.show()

if __name__ == '__main__': 
  #testOneTriangle()
  #testOneTriangle2()
  #testTwoTriangles()
  testRandomTriangles()



