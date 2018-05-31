#! /usr/bin/env python2.7
## -*- encoding: utf-8 -*-



## Usage
from sys import argv
if len(argv) < 4:
	from sys import stderr, exit
	stderr.write("Usage: {} $job[=$value] $npts ${{input}}.log [${{output}}]\n".format(argv[0]))
	exit(1)

from physics import A_to_a0
from orbkit import grid, core, read, output, extras
from numpy import amax, amin



## Input
from cclib.parser import ccopen
qc = read.convert_cclib(ccopen(argv[3]).parse(), all_mo=True)



## Get grid parameters and initialize grid
over_s = 7   # to be tuned

## Try to treat size parameter as a number
try:
	par = int(argv[2])
	## Two regimes: positive or zero, and negative
	## Positive or zero: the number of points is given and the spacing is deduced
	if 0 < par:
		grid.N_ = [par]*3

	elif par == 0:
		grid.N_ = [80]*3

	## Negative: the spacing is given and the number of points is deduced
	else:
		grid.delta_ = [1.0/(-par)]*3

	del par

## Didn't work - parameter is a keyword
except ValueError:
	grid.delta_ = [1.0/3.0  if argv[2] == "Coarse" else \
	               1.0/6.0  if argv[2] == "Medium" else \
	               1.0/12.0 if argv[2] == "Fine" else None]*3

grid.max_ = list(amax(qc.geo_spec.T, axis=1) + over_s)
grid.min_ = list(amin(qc.geo_spec.T, axis=1) - over_s)

grid.init()



## Get job type and parameters
#from orbkit import options
#options.numproc = 4

val = argv[1].split("=")
job = val[0]

if job == "MO":
	## Get list of orbitals
	mos = val[1].split(",")

	def func(data):
		return core.rho_compute(data, calc_mo=True, numproc=4,)

elif job == "FDensity":
	try:
		## Get keyword
		opt = val[1]
	except IndexError:
		## Default is SCF
		## (Useles for now, but keep it for later)
		opt = "SCF"
	mos = ["1:homo+1"]

	def func(data):
		return core.rho_compute(data, numproc=4)

elif job == "Potential":
	opt = val[1]



## Main calculation
qc.mo_spec = read.mo_select(qc.mo_spec, mos)["mo_spec"]
out = func(qc)



## Output
#try
	## Try outputting to file
#	output.main_output(out, qc.geo_info, qc.geo_spec, outputname=argv[4], otype="cb")

#except IndexError:
## No filename supplied - it's a visualisation
x, y, z = grid.x, grid.y, grid.z

## Holdover from the example code from which this was originally derived
## Set to True once mayavi works
mayavi_yes = True

if mayavi_yes:
	## If selected, use mayavi2 to make isosurface plot
	maya = False
	try:
		from enthought.mayavi import mlab
		maya = True
	except Exception:
		print("import enthought.mayavi failed -- trying mayavi")
	try:
		from mayavi import mlab
		maya = True
	except Exception:
		print("import mayavi failed")

	if maya:
		if len(qc.geo_spec) > 1:
			## Calculate best fitting plane via PCA
			from numpy import cov, mean, linalg
			from math import sqrt, atan2
			normal = linalg.eig(cov((qc.geo_spec - mean(qc.geo_spec, axis=0)).T))[1][:,-1]
			r_p = sqrt(normal[0]**2 + normal[1]**2)
			r = sqrt(r_p**2 + normal[2]**2)
			a, e = atan2(normal[0], normal[1]), atan2(r, r_p)
			mlab.view(azimuth=a, elevation=e, distance=min(grid.max_))

		mlab.figure(bgcolor=(1,1,1))
		for i, series in enumerate(out):

			mlab.contour3d(x, y, z, series, contours=[ 0.05 ], color=(0.4, 0, 0.235))
			mlab.contour3d(x, y, z, series, contours=[-0.05 ], color=(0.95, 0.95, 0.95))

			from numpy import stack
			mlab.points3d(*qc.geo_spec.T, extent=stack((grid.min_, grid.max_), axis=1).flatten(), color=(0,0,1))
			for p, t in stack((qc.geo_spec, qc.geo_info), axis=1):
				mlab.text3d(p[0], p[1], p[2], t[0], color=(0,0,0))

			mlab.show()

		#	mlab.savefig("./{}-{}.png".format(argv[4], i))
			mlab.clf()
