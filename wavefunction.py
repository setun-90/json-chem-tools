#! /usr/bin/env python2.7

## Structures for calculating atomic and molecular orbital wavefunctions.

from math import factorial as fac, fabs
from cmath import pi, sqrt, exp, sin, cos

## o_i(N) = S_s(c_si x_s(N))(M)
## x_s === G_s(r; a_s, R_s) = f(r) e^(-a_s|r - R_s|^2)
## f(r) = 

## Z === a_s : primitive_exponents|contraction_coefficients ?
## c_si      : (alpha|beta)_MO_coefficients ?

def binomial(n, m):
	try:
		return fac(n)//(fac(m)*fac(n-m))
	except ValueError:
		return 0

## GTO normalization constant
def N(a, lx, ly, lz):
	return sqrt((2**(2*(lx + ly + lz))*fac(lx)*fac(ly)*fac(lz)*pow(a, l + 1.5))
	           /(fac(2*lx)*fac(2*ly)*fac(2*lz)*pow(pi, 1.5)))

## Cartesian GTO angular part
def gc(lx, ly, lz, a, x, y, z):
	return N(a, lx, ly, lz) * x**lx*y**ly*z**lz * exp(-a*(x**2 + y**2 + z**2))

## Calculation of coefficient for transforming Cartesian GTO to Spherical GTO
## This can be memoized for an easy performance gain
def c(m, lx, ly, lz):
#	memo = {}
#	try:
#		return memo[(m,lx,ly,z)]
#	except KeyError:
	ma = abs(m)
	l = lx + ly + lz
	j = (lx + ly - ma)/2
#	memo[(m,lx,ly,lz)] = ...
	return sqrt(float((fac(2*lx)*fac(2*ly)*fac(2*lz)*fac(l)*fac(l - ma)) / (fac(2*l)*fac(lx)*fac(ly)*fac(lz)*(l + ma)))) * 1.0/(2**l*fac(l))  \
	     * sum(binomial(l, i)*binomial(i, j)*(-1)**i*fac(2*l - 2*i)/fac(l - ma - 2*i) for i in range(l - ma/2)) \
	     * sum(binomial(j, k)*binomial(ma, lx - 2*k) for k in range(j)) \
	     * (-1)**((ma - lx + 2*k)/2)

def v(m, lx, ly, lz, a, x, y, z):
	return c(m,lx,ly,lz)*gc(lx,ly,lz,a,x,y,z)

## Cartesian-to-Spherical GTO map
## c.f. H. B. Schlegel and M. J. Frisch, Int. J. of Quantum Chemistry, 54, 83-87 (1995)
vc = {
	(0, 0): lambda a, x, y, z: v( 0,0,0,0,a,x,y,z),

	(1,-1): lambda a, x, y, z: (v(1,1,0,0,a,x,y,z) - j*v(1,0,1,0,a,x,y,z))/sqrt(2.0),
	(1, 0): lambda a, x, y, z: v( 0,0,0,1,a,x,y,z),
	(1, 1): lambda a, x, y, z: (v(1,1,0,0,a,x,y,z) + j*v(1,0,1,0,a,x,y,z))/sqrt(2.0),

	(2,-2): lambda a, x, y, z: sqrt(3.0/8.0)*(v(2,2,0,0,a,x,y,z) - v(2,0,2,0,a,x,y,z)) - j*v(2,1,1,0,a,x,y,z)/sqrt(2.0),
	(2,-1): lambda a, x, y, z: (v(1,1,0,1,a,x,y,z) - j*v(1,0,1,1,a,x,y,z))/sqrt(2.0),
	(2, 0): lambda a, x, y, z: v(0,0,0,2,a,x,y,z) - (v(0,2,0,0,a,x,y,z) + v(0,0,2,0,a,x,y,z))/2.0,
	(2, 1): lambda a, x, y, z: (v(1,1,0,1,a,x,y,z) + j*v(1,0,1,1,a,x,y,z))/sqrt(2.0),
	(2, 2): lambda a, x, y, z: sqrt(3.0/8.0)*(v(2,2,0,0,a,x,y,z) - v(2,0,2,0,a,x,y,z)) + j*v(2,1,1,0,a,x,y,z)/sqrt(2.0),

	(3,-3): lambda a, x, y, z: sqrt(5.0)/4.0*(v(3,3,0,1,a,x,y,z) + j*v(3,0,3,0,a,x,y,z)) - 3.0/4.0*(v(3,1,2,0,a,x,y,z) + j*v(3,2,1,0,a,x,y,z)),
	(3,-2): lambda a, x, y, z: sqrt(3.0/8.0)*(v(2,2,0,1,a,x,y,z) - v(2,0,2,1,a,x,y,z)) - j*v(2,1,1,1,a,x,y,z)/sqrt(2.0), 
	(3,-1): lambda a, x, y, z: sqrt(3.0/5.0)*(v(1,1,0,2,a,x,y,z) - j*v(1,0,1,2,a,x,y,z)) - sqrt(3.0)*(v(1,3,0,0,a,x,y,z) - j*v(1,0,3,0,a,x,y,z))/4.0 - sqrt(3.0)*(v(1,1,2,0,a,x,y,z) - j*v(1,2,1,0,a,x,y,z))/(4.0*sqrt(5.0)), 
	(3, 0): lambda a, x, y, z: v(0,0,0,3,a,x,y,z) - 3.0*(v(0,2,0,1,a,x,y,z) + v(0,0,2,1,a,x,y,z))/(2.0*sqrt(5.0)),
	(3, 1): lambda a, x, y, z: sqrt(3.0/5.0)*(v(1,1,0,2,a,x,y,z) + j*v(1,0,1,2,a,x,y,z)) - sqrt(3.0)*(v(1,3,0,0,a,x,y,z) + j*v(1,0,3,0,a,x,y,z))/4.0 - sqrt(3.0)*(v(1,1,2,0,a,x,y,z) + j*v(1,2,1,0,a,x,y,z))/(4.0*sqrt(5.0)),
	(3, 2): lambda a, x, y, z: sqrt(3.0/8.0)*(v(2,2,0,1,a,x,y,z) - v(2,0,2,1,a,x,y,z)) + j*v(2,1,1,1,a,x,y,z)/sqrt(2.0),
	(3, 3): lambda a, x, y, z: sqrt(5.0)/4.0*(v(3,3,0,1,a,x,y,z) - j*v(3,0,3,0,a,x,y,z)) - 3.0/4.0*(v(3,1,2,0,a,x,y,z) - j*v(3,2,1,0,a,x,y,z)),

	(4,-3): lambda a, x, y, z: sqrt(5.0)/4.0*(v(3,3,0,1,a,x,y,z) + j*v(3,0,3,1,a,x,y,z)) - 3.0/4.0*(v(3,1,2,1,a,x,y,z) + j*v(3,2,1,1,a,x,y,z)),
	(4, 0): lambda a, x, y, z: v(0,0,0,4,a,x,y,z) + 3.0*(v(0,4,0,0,a,x,y,z) + v(0,0,4,0,a,x,y,z))/8.0 - 3.0*sqrt(3.0)*(v(0,2,0,2,a,x,y,z) + v(0,0,2,2,a,x,y,z) - 0.25*v(0,2,2,0,a,x,y,z))/sqrt(35.0),
	(4, 3): lambda a, x, y, z: sqrt(5.0)/4.0*(v(3,3,0,1,a,x,y,z) - j*v(3,0,3,1,a,x,y,z)) - 3.0/4.0*(v(3,1,2,1,a,x,y,z) - j*v(3,2,1,1,a,x,y,z)) #,

#	(5, 0): lambda a, x, y, z: 
}

## Spherical GTO angular part
def gs(l, m, a, x, y, z):
	try:
		return gc[(l,m)](a,x,y,z)
	except KeyError:
		raise ValueError("The quantum numbers describe an "
		               + ("impossible (l < |m|)" if l < abs(m) else "unimplemented")
		               + " atomic orbital")
