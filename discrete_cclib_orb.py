#! /usr/bin/env python2.7
## -*- encoding: utf-8 -*-



## Usage
from sys import argv
if len(argv) < 5:
	from sys import stderr, exit
	stderr.write("Usage: {} $action $job[=$value] $npts ${{input}}.log [${{output}}]\n".format(argv[0]))
	exit(1)

from physics import A_to_a0
from orbkit import grid, core, read, output, extras
from numpy import amax, amin



## Input
from cclib.parser import ccopen
qc = read.convert_cclib(ccopen(argv[4]).parse(), all_mo=True)



## Get grid parameters and initialize grid

## Oversizing in Bohr units
## ORBKIT uses 5 by default, tune this as required
#over_s = 7

## Spacing
grid.delta_ = [1.0/10.0]*3

#grid.max_ = list(amax(qc.geo_spec.T, axis=1) + over_s)
#grid.min_ = list(amin(qc.geo_spec.T, axis=1) - over_s)
grid.max_ = [ 20]*3
grid.min_ = [-20]*3

grid.init()



## Get job type and parameters
val = argv[2].split("=")
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
if argv[1] == "save":
	## Save to file
	try:
		output.main_output(out, qc.geo_info, qc.geo_spec, outputname=argv[4], otype="cb")

	except IndexError:
		from sys import stderr, exit
		stderr.write("No filename specified.\n")
		exit(1)

elif argv[1] == "viz":
	## Vizualize with MayaVi
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

			from numpy import stack, column_stack
			extents = column_stack((grid.min_, grid.max_)).flatten()
			mlab.contour3d(series, contours=[ 0.05 ], extent=extents, color=(0.4, 0, 0.235))
			mlab.contour3d(series, contours=[-0.05 ], extent=extents, color=(0.95, 0.95, 0.95))

			from IPython import embed as shell

			mlab.points3d(*(qc.geo_spec.T/A_to_a0), color=(0,0,0), scale_factor=1)
		#	for p, t in stack((qc.geo_spec, qc.geo_info), axis=1):
		#		mlab.text3d(p[0], p[1], p[2], t[0], color=(0,0,0))

			mlab.show()

			shell()

			mlab.savefig("./{}-{}.png".format(argv[4], i))
			mlab.clf()
