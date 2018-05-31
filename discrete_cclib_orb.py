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

grid.max_ = map(lambda a: max(a) + over_s, qc.geo_spec.T)
grid.min_ = map(lambda a: min(a) - over_s, qc.geo_spec.T)

grid.init()



## Get job type and parameters
#from orbkit import options
#options.numproc = 4

val = argv[1].split("=")
job = val[0]

if job == "MO":
	## Get list of orbitals
	mos = val[1].split(",")

	qc.mo_spec = read.mo_select(qc.mo_spec, mos)["mo_spec"]

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
	qc.mo_spec = read.mo_select(qc.mo_spec, mos)["mo_spec"]	

	def func(data):
		return core.rho_compute(data, numproc=4)

elif job == "Potential":
	opt = val[1]



## Main calculation
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
		## Calculate best fitting plane
		if len(qc.geo_spec) > 1:
			## PCA
			from numpy import cov, mean, linalg
			from math import sqrt, atan2
			eival, eivec = linalg.eig(cov((qc.geo_spec - mean(qc.geo_spec, axis=0)).T))
			normal, point = eivec[:,-1], mean(qc.geo_spec, axis=0)
			r_p = sqrt(normal[0]**2 + normal[1]**2)
			r = sqrt(r_p**2 + normal[2]**2)
			a, e = atan2(normal[0], normal[1]), atan2(r, r_p)
			mlab.view(azimuth=a, elevation=e)

		mlab.figure(bgcolor=(1,1,1))
		for i, series in enumerate(out):
			P_x, P_y, P_z = qc.geo_spec.T
			mlab.points3d(P_x, P_y, P_z, color=(0,0,1), mode="sphere", scale_factor=1)

			src = mlab.pipeline.scalar_field(series)
			mlab.pipeline.iso_surface(src, contours=[ 0.05, ], opacity=1, color=(0.4, 0, 0.235))
			mlab.pipeline.iso_surface(src, contours=[-0.05, ], opacity=1, color=(0.95, 0.95, 0.95))
			
			mlab.savefig("./{}-{}.png".format(argv[4], i))
			mlab.clf()

else:
	## Use matplotlib to show cuts of the molecular orbitals
	import matplotlib.pyplot as plt
	import numpy as np

	## Select cuts
	xd = out[0][grid.N_[0]/2-1,:,:]
	yd = out[0][:,grid.N_[1]/2-1,:]
	zd = out[0][:,:,grid.N_[2]/2-1]

	## Plot cuts
	f, (pic1, pic2, pic3) = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(6,14))
	pic1.contour(z, y, xd, 50, linewidths=0.5, colors='k')
	pic1.contourf(z, y, xd, 50, cmap=plt.cm.rainbow, vmax=abs(xd).max(), vmin=-abs(xd).max())
	pic1.set_xlabel('z')
	pic1.set_ylabel('y')

	pic2.contour(z, x, yd, 50, linewidths=0.5, colors='k')
	pic2.contourf(z, x, yd, 50, cmap=plt.cm.rainbow, vmax=abs(yd).max(), vmin=-abs(yd).max())
	pic2.set_xlabel('z')
	pic2.set_ylabel('x')

	pic3.contour(y, x, zd, 50, linewidths=0.5, colors='k')
	pic3.contourf(y, x, zd, 50, cmap=plt.cm.rainbow, vmax=abs(zd).max(), vmin=-abs(zd).max())
	pic3.set_xlabel('y')
	pic3.set_ylabel('x')

	## Following options applied for all subplots as they share x- and y-axis
	pic1.xaxis.set_ticks(np.arange(-5,6,5))
	pic1.yaxis.set_ticks(np.arange(-5,6,5))
	pic1.set_aspect('equal')

	## Plot
	f.subplots_adjust(left=0.15,bottom=0.05,top=0.95,right=0.95)
	f.show()

	raw_input("Press Enter to continue")
