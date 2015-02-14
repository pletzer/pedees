#!/usr/bin/env python

import numpy
import math

class Cg:

	def __init__(self, mat, b):
		self.mat = mat
		self.b = b

	def solve(self, precond, x0, numIters, tol, verbose = True):

		x = x0[:]
		r = self.b - self.matrixDotVector(x)
		z = r / precond
		p = z[:]

		k = 0
		err = self.norm(r)
		while k < numIters:
			ap = self.matrixDotVector(p)
			alpha = numpy.dot(r, z) / numpy.dot(p,  ap)
			x += alpha * p
			rOld = r[:]
			r -= alpha * ap
			err = self.norm(r)

			if err < tol:
				break

			z = r / precond
			beta = numpy.dot(z, r)/numpy.dot(z, rOld)
			p = z + beta*p
			k += 1

		return x

	def matrixDotVector(self, vec):
		res = numpy.zeros( vec.shape, numpy.float64 )
		for ij, matVal in self.mat.items():
			i, j = ij
			res[i] += matVal * vec[j]
		return res

	def norm(self, v):
		return math.sqrt( numpy.dot(v, v) )

###############################################################################
def test1():
	mat = {(0,0): 1.0}
	b = numpy.array([1.0])
	s = Cg(mat, b)
	p = numpy.array([1.0])
	x = s.solve(precond = p, x0 = numpy.array([0.0]), numIters=10, tol=1.e-10, verbose=True)
	print 'x = ', x

def test3():
	mat = {(0,0): 1.0, (1,1): 1.0, (2,2): 1.0}
	b = numpy.array([1.0, 2.0, 3.0])
	s = Cg(mat, b)
	p = numpy.array([1.0, 1.0, 1.0])
	x = s.solve(precond = p, x0 = numpy.array([0.0, 0.0, 0.0]), 
		numIters=10, tol=1.e-10, verbose=True)
	print 'x = ', x

def test3b():
	mat = {(0,0): 1.0, (1,1): 2.0, (2,2): 3.0}
	b = numpy.array([1.0, 2.0, 3.0])
	s = Cg(mat, b)
	p = numpy.array([1.0, 2.0, 3.0])
	x = s.solve(precond = p, x0 = numpy.array([0.0, 0.0, 0.0]), 
		numIters=10, tol=1.e-10, verbose=True)
	print 'x = ', x


if __name__ == '__main__':
	test1()
	test3()
	test3b()