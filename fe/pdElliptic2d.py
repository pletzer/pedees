#!/usr/bin/env python

class Elliptic2d:

	def __init__(self, fFunc, gFunc, sFunc):

		self.fFunc = fFunc
		self.gFunc = gFunc
		self.sFunc = sFunc

		self.mat = {}
		self.sourceVec = []

	def assemble(self, points, triangles):

		n = len(points)
		self.sourceVec = numpy.zeros( (n,), numpy.float64 )

		for iabc in triangles:

			ia, ib, ic = iabc
			pa, pb, pc = points[ia], points[ib], points[ic]

			# centroid 
			pMid = (pa + pb + pc)/3.0
			fxx = fyy = self.fFunc(pMid)

			ga = self.gFunc(pa)
			gb = self.gFunc(pb)
			gc = self.gFunc(pc)

			sa = self.sFunc(pa)
			sb = self.sFunc(pb)
			sc = self.sFunc(pc)

			xcb = pc[0] - pb[0]
			ycb = pc[1] - pb[1]
			xac = pa[0] - pc[0]
			yac = pa[1] - pc[1]
			xba = pb[0] - pa[0]
			yba = pb[1] - pa[1]

			area = -xba*yac + yba*xac

			fOverA = 0.5*(fxx + fyy)/area

			faa = fOverA * (ycb*ycb + xcb*xcb) \
			      + (ga/ 20. + gb/ 60. + gc/ 60.)*area

			fab = fOverA * (ycb*yac + xcb*xac) \
			      + (ga/ 60. + gb/ 60. + gc/120.)*area

			fac = fOverA * (ycb*yba + xcb*xba) \
			      + (ga/ 60. + gb/120. + gc/ 60.)*area

			fbb = fOverA * (yac*yac + xac*xac) \
			      + (ga/ 60. + gb/ 20. + gc/ 60.)*area

			fbc = fOverA * (yac*yba + xac*xba) \
			      + (ga/120. + gb/ 60. + gc/ 60.)*area

			fcc = fOverA * (yba*yba + xba*xba) \
			      + (ga/ 60. + gb/ 60. + gc/ 20.)*area

			self.mat[ia, ia] = self.mat.get((ia, ia), 0.0) + faa
			self.mat[ia, ib] = self.mat.get((ia, ib), 0.0) + fab
			self.mat[ia, ic] = self.mat.get((ia, ic), 0.0) + fac
			self.mat[ib, ib] = self.mat.get((ib, ib), 0.0) + fbb
			self.mat[ib, ic] = self.mat.get((ib, ic), 0.0) + fbc
			self.mat[ic, ic] = self.mat.get((ic, ic), 0.0) + fcc

			# make sure matrix is Hermitian
			self.mat[ib, ia] = self.mat[ia, ib]
			self.mat[ic, ia] = self.mat[ia, ic]
			self.mat[ic, ib] = self.mat[ib, ic]

			self.sourceVec[ia] += area*(sa/12.0 + sb/24.0 + sc/24.0)
			self.sourceVec[ib] += area*(sa/24.0 + sb/12.0 + sc/24.0)
			self.sourceVec[ic] += area*(sa/24.0 + sb/24.0 + sc/12.0)

	def applyBoundaryConditions(self, edges, aVals, bVals, cVals):
		pass

	def getStiffnessMatrix(self):
		return self.mat

	def getSourceVector(self):
		return self.sourceVec