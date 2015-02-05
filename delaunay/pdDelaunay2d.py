#!/usr/bin/env python

import numpy
import math

class Delaunay2d:

  # threshold for declaring a triangle non-degenerate
  self.MIN_AREA = 1.23456789e-14
  self.EPS = 1.23456789e-14

	def __init__(self, points, edges, holes):
    """
    Constructor 
    @param points list of points
    @param edges list of boundary edges, anticlockwise for the out edges, clockwise
                  for the inner edges
    @param holes a coordinate point inside each hole
    """

    # list of vertices
		self.xyPoints = points[:]

    # list of integer 3-tuples
    self.triangles = []

    # list of boundary edges, so far
    self.boundaryEdges = []

    # map from edge 2-tuples to triangle 3-tuples
    edges2Triangles = {}

    # not currently used
    self.edges = edges
    self.holes = holes

	def triangulate(self):
    pass

  def getTriangles(self):
    pass

  def getEdges(self):
    pass

  def _getArea(self, ip0, ip1, ip2):
    """
    Compute the area spanned by the traingle ip0, ip1, and ip2
    @param ip0 index of first vertex
    @param ip1 index of second vertex
    @param ip2 index of third vertex
    """
    d1 = self.points[ip1] - self.points[ip0]
    d2 = self.points[ip2] - self.points[ip0]
    return 0.5 * 0.5*(d1[0]*d2[1] - d1[1]*d2[0])

  def _isEdgeVisible(self, ip, edge):
    """
    Return true if edge is visible by a point. An edge is visible 
    iff the point lies to its right. 
    @param p point index
    @param edge (2 point indices)
    @return True if visible    
    """
    area = self._getArea(ip, edge[0], edge[1])
    if area > self.MIN_AREA:
      return True
    return False

  def _addPoint(self, ip):
    """
    Add point
    @param ip point index
    """
    # For each point we determine which 
    # boundary edges are visible to the 
    # new point. Create a new triangle between point and 
    # visible edge. Collect the new edges and update the 
    # list of triangles, the edge to triangle connectivity, 
    # and the list of boundary edges.

    edgeIndicesToRemove = []
    newBoundaryEdges = set()

    for ibedge in range(0, len(self.boundaryEdges):
      # iterate over the convex hull (boundary) edges 

      boundEdge = self.boundaryEdges[ibedge];

      # an edge is visible from the new point "ipoint" iff
      # the point is to the right of the boundary, which is
      # assumed to go counterclockwise.
      if self.isEdgeVisible(ipoint, boundEdge):
        # add new triangle
        newT = [boundEdge[0], ipoint, boundEdge[1]]
        self.triangles.append(newT)

        # tag this boundary edge for removal
        edgeIndicesToRemove.append(ibedge)

        # update the edge to triangles map
        it = self.edge2Triangles[boundEdge]
        itriangle = len(self.triangles) - 1
        it.append(itriangle)

        # cache the new boundary edges
        e1 = (newT[0], newT[1])
        e2 = (newT[1], newT[2])
        newBoundaryEdges.append(e1)
        newBoundaryEdges.append(e2)

        # add the edge -> connectivity entries
        lastTriangleIndex = len(self.triangles) - 1
        newTriangleVectPair = [lastTriangleIndex]
        self.edge2Triangles[e1] = newTriangleVectPair
        self.edge2Triangles[e2] = newTriangleVectPair

    # remove the old boundary edges
    for i in range(len(edgeIndicesToRemove) - 1), -1, -1):
      ie = edgeIndicesToRemove[i];
      del self.boundaryEdges[self.boundaryEdges[ie]

    # add the new boundary edges (i, j), but only if (j, i) is 
    # not in the list. If both (i, j) and (j, i) are in the list
    # then this is not a boundary edge. 
    for edge in newBoundaryEdges:
      complEdge = (edge[1], edge[0])
      if newBoundaryEdges.has_key(complEdge):
        # complEdge was not found so add edge
        self.boundaryEdges.append(edge)
      else:
        # complEdge was found
        # merge the edge to triangle index connectivity map
        e2t = self.edge2Triangles.get(edge, None)
        complE2t = self.edge2Triangles.get(complEdge, None)
        if e2t is not None:
          e2t.append(complE2t[0])
          del self.edge2Triangles[complEdge]

def _flipEdge(self, bedge):

  tris = self.edge2Triangles.get(bedge)

  # must have two adjacent triangles
  if tris == None or len(tris) < 2: 
    return False

  t1, t2 = tris[0], tris[1]
  triangle1, triangle2 = self.triangles[t1], self.triangles[t2]
  # opposite nodes i1 and i2
  i1, i2 = -1, -1
  for i in range(3):
    found1, found2 = False, False
    for j in range(2):
      found1 |= (bedge[j] == triangle1[i])
      found2 |= (bedge[j] == triangle2[i])
    if not found1: i1 = i
    if not found2: i2 = i

    # get the coordinates of the quadrilateral
    #     b
    #    /^\
    #   / | \
    #  /  |  \
    # c 1 | 2 d
    #  \  |  /
    #   \ | /
    #    \|/ 
    #     a
    pa, pb = bedge[0], bedge[1]
    pc, pd = triangle1[i1], triangle2[i2]
    xya, xyb, xyc, xyd = self.xyPoints[pa], self.xyPoints[pb], self.xyPoints[pc], self.xyPoints[pd]

    # the 2 angles opposite to bedge
    crossProd1 = 2.0*self.getArea(pc, pa, pb)
    crossProd2 = 2.0*self.getArea(pd, pb, pa)
    dotProd1 = (xya[0]-xyc[0])*(xyb[0]-xyc[0]) + (xya[1]-xyc[1])*(xyb[1]-xyc[1])
    dotProd2 = (xyb[0]-xyd[0])*(xya[0]-xyd[0]) + (xyb[1]-xyd[1])*(xya[1]-xyd[1])
    angle1 = abs(math.atan2(crossProd1, dotProd1))
    angle2 = abs(math.atan2(crossProd2, dotProd2))
   
    # Delaunay test
    if angle1 + angle2 > math.pi*(1.0 + self.EPS):
      # flip
      #     b
      #    / \
      #   /   \
      #  /  1  \
      # c-----> d
      #  \  2  /
      #   \   /
      #    \ / 
      #     a
      newTriangle1 = [pc, pd, pb]
      newTriangle2 = [pd, pc, pa]
      self._makeCounterClockwise(newTriangle1)
      self._makeCounterClockwise(newTriangle2)

      # reset the triangle list
      self.triangles[t1] = newTriangle1
      self.triangles[t2] = newTriangle2

      # create new entry in edge -> triangle map
      newEdge = (pc, pd)
      t12 = [t1, t2]
      self.edge2Triangles[newEdge] = t12

      # update neighboring entries
      otherEdge = (pb, pd)
      complOtherEdge = (pd, pb)
      if self.edge2Triangles.has_key(otherEdge):
        index = 0

      otherEdge[0] = pb; otherEdge[1] = pd
      complOtherEdge[0] = pd; complOtherEdge[1] = pb
      otherEdgeIt = self.edge2Triangles.get(otherEdge, None)
      complOtherEdgeIt = self.edge2Triangles.get(complOtherEdge, None)
      if otherEdgeIt != None:
        index = 0
        if len(otherEdgeIt) > 1 and otherEdgeIt[1] == t2:
          index = 1
        otherEdgeIt[index] = t1
      elif complOtherEdgeIt != None:
        index = 0
        if len(complOtherEdgeIt) > 1 and complOtherEdgeIt[1] == t2:
          index = 1
        complOtherEdgeIt[index] = t1

      otherEdge[0] = pc; otherEdge[1] = pa
      complOtherEdge[0] = pa; complOtherEdge[1] = pc
      otherEdgeIt = self.edge2Triangles.get(otherEdge, None)
      complOtherEdgeIt = self.edge2Triangles.get(complOtherEdge, None)
      if otherEdgeIt != None:
        index = 0
        if len(otherEdgeIt) > 1 and otherEdgeIt[1] == t1:
          index = 1
        otherEdgeIt[index] = t2
      elif complOtherEdgeIt != None:
        index = 0
        if (len(complOtherEdgeIt) and complOtherEdgeIt[1] == t1:
          index = 1
        complOtherEdgeIt[index] = t2

      # remove entry
      del self.edge2Triangles[bedge]

      return True
    
    return False

######################################################

def test():
  

if __name__ == '__main__': test()