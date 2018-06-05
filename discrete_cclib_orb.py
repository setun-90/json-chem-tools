#! /usr/bin/env python2.7
## -*- encoding: utf-8 -*-



## Usage
from sys import argv
if len(argv) < 5:
	from sys import stderr, exit
	stderr.write("Usage: {} $action $job[=$value] $npts ${{input}}.log ${{output}}\n".format(argv[0]))
	exit(1)

from physics import A_to_a0
from orbkit import grid, core, output, extras
## Patched version of orbkit.read
import read
from numpy import amax, amin



## Input
if ".json" in argv[4]:
	from json import load
	qc = read.convert_json(load(open("fchk_log_files/cid241/OPT.json", "r")), all_mo=True)
else:
	from cclib.parser import ccopen
	qc = read.convert_cclib(ccopen(argv[4]).parse(), all_mo=True)



## Get grid parameters and initialize grid

## Oversizing in Bohr radii
## ORBKIT uses 5 by default, tune this as required
over_s = 7

## Spacing/Number of points
par = int(argv[3])
if par > 0:
	grid.N_ = [par]*3
elif par == 0:
	grid.N_ = [80]*3
else:
	grid.delta_ = [1.0/(-par)]*3

grid.max_ = amax(qc.geo_spec.T, axis=1) + over_s
grid.min_ = amin(qc.geo_spec.T, axis=1) - over_s

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
		for i, series in enumerate(out):
			output.main_output(series, qc.geo_info, qc.geo_spec, outputname="{}-{}".format(argv[5], i), otype="cb")

	except IndexError:
		from sys import stderr, exit
		stderr.write("No output filename specified.\n")
		exit(1)

elif argv[1] == "viz":
	## Vizualize with MayaVi
	maya = False
	try:
		from enthought.mayavi import mlab
		maya = True
	except ImportError as E:
		print("import enthought.mayavi failed -- trying mayavi")
	try:
		from mayavi import mlab
		maya = True
	except ImportError:
		print("import mayavi failed")

	if maya:
		mlab.figure(bgcolor=(1,1,1))
		if len(qc.geo_spec) > 1:
			## Calculate best fitting plane via PCA
			from numpy import cov, mean, linalg
			from math import sqrt, atan2
			normal = linalg.eig(cov((qc.geo_spec - mean(qc.geo_spec, axis=0)).T))[1][:,-1]
			r_p = sqrt(normal[0]**2 + normal[1]**2)
			r = sqrt(r_p**2 + normal[2]**2)
			a, e = atan2(normal[0], normal[1]), atan2(r, r_p)
			mlab.view(azimuth=a, elevation=e, distance=min(grid.max_))

		for i, series in enumerate(out):

			from numpy import mgrid
			X, Y, Z = mgrid[grid.min_[0]:grid.max_[0]:1j*len(grid.x),
			                grid.min_[1]:grid.max_[1]:1j*len(grid.y),
			                grid.min_[2]:grid.max_[2]:1j*len(grid.z)]

			data = mlab.pipeline.scalar_field(X, Y, Z, series)

			op = mlab.pipeline.iso_surface(data, contours=[ 0.05 ], color=(0.4, 0, 0.235))
			op.actor.property.interpolation = "phong"
			op.actor.property.specular = 0.1
			op.actor.property.specular_power = 5
			on = mlab.pipeline.iso_surface(data, contours=[-0.05 ], color=(0.95, 0.95, 0.95))
			on.actor.property.interpolation = "phong"
			on.actor.property.specular = 0.1
			on.actor.property.specular_power = 5


			tab = []
		#	mlab.points3d(*qc.geo_spec.T, color=(0,0,0), scale_factor=1)
			with open("Atoms.csv", "r") as f:
				f.next() # skip header line
				tab = [line.split() for line in f]
			tran = dict([(line[1], line[2:]) for line in tab])

			from numpy import stack
			for P, t in stack((qc.geo_spec, qc.geo_info), axis=1):
				attr = tran[t[0]]
				color = tuple(map(lambda x: float(x)/255.0, attr[1:4]))
				scale_factor = float(attr[4])*0.4
				mlab.points3d([float(P[0])], [float(P[1])], [float(P[2])], color=color, scale_factor=scale_factor, resolution=15)
		#		mlab.text3d(p[0], p[1], p[2], t[0], color=(0,0,0))

			mlab.show()

			from IPython import embed as shell
		#	shell()

			mlab.savefig("./{}-{}.png".format(argv[5], i))
			mlab.clf()
