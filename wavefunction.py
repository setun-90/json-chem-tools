## -*- encoding: utf-8 -*-



## Structures for calculating atomic and molecular orbital wavefunctions.
from math import factorial as fac, fabs
from cmath import pi, sqrt, exp, sin, cos



## Memoization
def memoize(f):
	memo = {}
	def aux(*x):
		try:
			return memo[x]
		except KeyError:
			memo[x] = f(*x)
			return memo[x]
	return aux



def binomial(n, m):
	try:
		return fac(n)//(fac(m)*fac(n-m))
	except ValueError:
		return 0

## GTO normalization constant
## Can be memoized because it doesn't depend on position
@memoize
def N(a, lx, ly, lz):
	l = lx + ly + lz
	return sqrt(float((2**(2*l)*fac(lx)*fac(ly)*fac(lz)*pow(a, l + 1.5))
	                 /(fac(2*lx)*fac(2*ly)*fac(2*lz)*pow(pi, 1.5))))

## Cartesian GTO angular part (normalized & contracted)
def gc(lx, ly, lz, bs, x, y, z):
	return x**lx*y**ly*z**lz * sum(cc[1] * N(cc[0], lx, ly, lz) * exp(-cc[0]*(x**2 + y**2 + z**2)) for cc in bs)

## Calculation of coefficient for transforming Cartesian GTO to Spherical GTO
## Can be memoized because it doesn't depend on position
@memoize
def c(m, lx, ly, lz):
	ma = abs(m)
	l = lx + ly + lz
	j = (lx + ly - ma)/2
	return sqrt(float((fac(2*lx)*fac(2*ly)*fac(2*lz)*fac(l)*fac(l - ma)) / (fac(2*l)*fac(lx)*fac(ly)*fac(lz)*fac(l + ma)))) * 1.0/(2**l*fac(l))  \
	     * sum(binomial(l, i)*binomial(i, j) * (-1)**i*fac(2*l - 2*i)/fac(l - ma - 2*i) for i in range((l - ma)//2 + 1)) \
	     * sum(binomial(j, k)*binomial(ma, lx - 2*k)*(-1)**((ma - lx + 2*k)/2) for k in range(j + 1))

def v(m, lx, ly, lz, bs, x, y, z):
	return c(m,lx,ly,lz) * gc(lx,ly,lz,bs,x,y,z)

## Cartesian-to-Spherical GTO map
## c.f. H. B. Schlegel and M. J. Frisch, Int. J. of Quantum Chemistry, 54, 83-87 (1995)
vc = {
	(0, 0): lambda bs, x, y, z: v(0,0,0,0,bs,x,y,z),

	(1, 0): lambda bs, x, y, z: v(0,0,0,1,bs,x,y,z),
	(1, 1): lambda bs, x, y, z: sqrt(2.0)*v(1,1,0,0,bs,x,y,z),

	(2, 0): lambda bs, x, y, z: v(0,0,0,2,bs,x,y,z) - (v(0,2,0,0,bs,x,y,z) + v(0,0,2,0,bs,x,y,z))/2.0,
	(2, 1): lambda bs, x, y, z: sqrt(2.0)*v(1,1,0,1,bs,x,y,z),
	(2, 2): lambda bs, x, y, z: 2.0*sqrt(3.0/8.0)*(v(2,2,0,0,bs,x,y,z) - v(2,0,2,0,bs,x,y,z)),

	(3, 0): lambda bs, x, y, z: v(0,0,0,3,bs,x,y,z) - 3.0*(v(0,2,0,1,bs,x,y,z) + v(0,0,2,1,bs,x,y,z))/(2.0*sqrt(5.0)),
	(3, 1): lambda bs, x, y, z: 2.0*(sqrt(3.0/5.0)*v(1,1,0,2,bs,x,y,z) - sqrt(3.0)*v(1,3,0,0,bs,x,y,z)/4.0 - sqrt(3.0)*v(1,1,2,0,bs,x,y,z)/(4.0*sqrt(5.0))),
	(3, 2): lambda bs, x, y, z: 2.0*sqrt(3.0/8.0)*(v(2,2,0,1,bs,x,y,z) - v(2,0,2,1,bs,x,y,z)),
	(3, 3): lambda bs, x, y, z: (sqrt(5.0)*v(3,3,0,1,bs,x,y,z) - 3.0*v(3,1,2,0,bs,x,y,z))/2.0,

	(4, 0): lambda bs, x, y, z: v(0,0,0,4,bs,x,y,z) + 3.0*(v(0,4,0,0,bs,x,y,z) + v(0,0,4,0,bs,x,y,z))/8.0 - 3.0*sqrt(3.0)*(v(0,2,0,2,bs,x,y,z) + v(0,0,2,2,bs,x,y,z) - 0.25*v(0,2,2,0,bs,x,y,z))/sqrt(35.0),
	(4, 1): lambda bs, x, y, z: 2.0*(sqrt(5.0/7.0)*v(1,1,0,3,bs,x,y,z) - 3.0*(sqrt(5.0)*v(1,3,0,1,bs,x,y,z) - v(1,1,2,1,bs,x,y,z))),
	(4, 2): lambda bs, x, y, z: 2.0*(sqrt(27.0/56.0)*(v(2,2,0,2,bs,x,y,z) - v(2,0,2,2,bs,x,y,z)) - sqrt(5.0/32.0)*(v(2,4,0,0,bs,x,y,z) - v(2,0,4,0,bs,x,y,z))),
	(4, 3): lambda bs, x, y, z: (sqrt(5.0)*v(3,3,0,1,bs,x,y,z) - 3.0*v(3,1,2,1,bs,x,y,z))/4.0,
	(4, 4): lambda bs, x, y, z: 2.0*(sqrt(35.0/128.0)*(v(4,4,0,0,bs,x,y,z) + v(4,0,4,0,bs,x,y,z)) - sqrt(27.0/32.0)*v(4,2,2,0,bs,x,y,z)) #,

#	(5, 0): lambda bs, x, y, z: (0,0,5) - ((2,0,3) + (0,2,3)) + ((4,0,1) + (0,4,1)) + (2,2,1),
#	(5, 1): lambda bs, x, y, z: 2.0*((1,0,4) - (3,0,2) - (1,2,2) + (5,0,0) + (1,4,0) + (3,2,0)),
#	(5, 2): lambda bs, x, y, z: 2.0*(((2,0,3) - (0,2,3)) - ((4,0,1) - (0,4,1))),
#	(5, 3): lambda bs, x, y, z: 2.0*((3,0,2) - (1,2,2) - ((5,0,0) -(1,4,0)) + (3,2,0)),
#	(5, 4): lambda bs, x, y, z: 2.0*(((4,0,1) + (0,4,1)) - (2,2,1)),
#	(5, 5): lambda bs, x, y, z: 2.0*((5,0,0) + (1,4,0) - (3,2,0)),
}

## Spherical GTO angular part
def gs(l, m, bs, x, y, z):
	try:
		return vc[l, m](bs, x, y, z).real
	except KeyError:
		raise ValueError("The quantum numbers describe an "
		               + ("impossible (l < |m|)" if l < abs(m) else "unimplemented")
		               + " atomic orbital")

## MO wavefunction
def psi_MO(basis_set, MO_coeffs, x, y, z):
	return [sum(C*gs(l, m, basis_set[1], x, y, z) for C, l in zip(arr, basis_set[0]) for m in range(l + 1)) for arr in MO_coeffs]

## Extraction of primitive exponents for each shell
## c.f. CHK-JSON-Shell.pdf
def primitives(basis_set):
	res = []
	for i, atom in enumerate(basis_set, 1):
		atomi = iter(zip(atom, atom[1:]))   # Necessary for using next()
		for ssh, n_ssh in atomi:
			if ssh[0] == "S" and n_ssh[0] == "P" and ssh[1][0][0] == n_ssh[1][0][0]:
				## It's an SP hybridization
				ssh[0] = "SP"
				for e, c in zip(ssh[1], n_ssh[1]):
					e += [c[1]]
				next(atomi, None)
			else:
				## It's a normal contraction
				for e in ssh[1]:
					e += [0]
			res += [(i, ssh)]
		## Last subshell
		for e in atom[-1][1]:
			e += [0]
		res += [(i, atom[-1])]

	return res
