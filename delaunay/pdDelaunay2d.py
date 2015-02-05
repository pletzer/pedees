#!/usr/bin/env python

TO DEBUG...

import numpy

class Delaunay2d:

  # threshold for declaring a triangle non-degenerate
  self.MIN_AREA = 1.23456789e-14

	def __init__(self, points, edges, holes):
    """
    Constructor 
    @param points list of points
    @param edges list of boundary edges, anticlockwise for the out edges, clockwise
                  for the inner edges
    @param holes a coordinate point inside each hole
    """
		self.points = points
    self.boundaryEdges = edges
    self.holes = holes
    self.triangles = []
    edges2Triangles = {}

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
        newT = [0, 0, 0]
        newT[0] = boundEdge[0]
        newT[1] = ipoint
        newT[2] = boundEdge[1]
        self.triangles.append(newT)

        # tag this boundary edge for removal
        edgeIndicesToRemove.append(ibedge)

        # update the edge to triangles map
        it = self.edge2Triangles[boundEdge]
        itriangle = self.triangles.size() - 1
        it.append(itriangle)

        # cache the new boundary edges
        e1(2) = [0, 0]
        e2(2) = [0, 0]
        e1[0] = newT[0]; e1[1] = newT[1]
        e2[0] = newT[1]; e2[1] = newT[2]
        newBoundaryEdges.append(e1);
        newBoundaryEdges.append(e2);

        # add the edge -> connectivity entries
        lastTriangleIndex = len(self.triangles) - 1
        newTriangleVectPair = [lastTriangleIndex]
        self.edge2Triangles[e1] = newTriangleVectPair
        self.edge2Triangles[e2] = newTriangleVectPair

    # remove the old boundary edges
    for i in range(len(edgeIndicesToRemove) - 1), -1):
      ie = edgeIndicesToRemove[i];
      del self.boundaryEdges[self.boundaryEdges[ie]

    # add the new boundary edges (i, j), but only if (j, i) is 
    # not in the list. If both (i, j) and (j, i) are in the list
    # then this is not a boundary edge. 
    for edge in newBoundaryEdges:
      complEdge = [-1, -1]
      complEdge[0] = edge[1]; complEdge[1] = edge[0]
      if newBoundaryEdges.has_key(complEdge):
        # complEdge was not found so add edge
        self.boundaryEdges.append(edge)
      else:
        # complEdge was found
        # merge the edge to triangle index connectivity map
        e2t = self.edge2Triangles.get(edge, None)
        complE2t = self.edge2Triangles.get(complEdge, None)
        if e2t != None:
          e2t.append(complE2t[0])
          del self.edge2Triangles[complEdge]

######################################################

def test():
  

if __name__ == '__main__': test()