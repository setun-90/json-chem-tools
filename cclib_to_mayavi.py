#! /usr/bin/env python2.7
## -*- encoding: utf-8 -*-



## Usage
from sys import argv
if len(argv) < 2:
	from sys import stderr, exit
	stderr.write("Usage: {} $npts ${{input}}".format(argv[0]))
	exit(1)

mayavi_yes = True

from orbkit import grid, core, read, output, extras



## Input
from cclib.parser import ccopen
qc = read.convert_cclib(ccopen(argv[2]).parse(), all_mo=True)



## Get grid parameters and initialize grid
over_s = 5   # to be tuned

## Try to treat size parameter as a number
try:
	p_npts = int(argv[1])
	## Two regimes: positive or zero, and negative
	## Positive or zero: the number of points is given and the spacing is deduced
	if 0 < p_npts:
		grid.N_ = [p_npts]*3

	elif p_npts == 0:
		grid.N_ = [80]*3

	## Negative: the spacing is given and the number of points is deduced
	else:
		## -1 is not implemented
		grid.delta_ = [2.0**(2+p_npts)/3.0   if -5 < p_npts < -1 else \
		               -p_npts*1e-3/A_to_a0  if p_npts <= -5 else None]*3

	del p_npts

## Didn't work - parameter is a keyword
except ValueError:
	grid.delta_ = [1.0/3.0  if argv[2] == "Coarse" else \
	               1.0/6.0  if argv[2] == "Medium" else \
	               1.0/12.0 if argv[2] == "Fine" else None]*3
#grid.delta_ = [1.0/6.0]*3
grid.max_ = list(map(lambda a: max(a) + over_s, qc.geo_spec))
grid.min_ = list(map(lambda a: min(a) - over_s, qc.geo_spec))
grid.init()



## Define the molecular orbital to be calculated
selected_MO = ['homo']
qc.mo_spec = read.mo_select(qc.mo_spec, selected_MO)["mo_spec"]



## Main calculation
mo_list = core.rho_compute(qc, calc_mo=True, numproc=4)



## Plot the results
x, y, z = grid.x, grid.y, grid.z

if mayavi_yes:
	## If selected, use mayavi2 to make isosurface plot
	maya = False
	try:
		from enthought.mayavi import mlab
		maya = True
	except Exception:
		pass
	try:
		from mayavi import mlab
		maya = True
	except Exception:
		pass

	if not maya:
		print("error importing mayavi")
	else:
		src = mlab.pipeline.scalar_field(mo_list[0])
		mlab.pipeline.iso_surface(src, contours=[0.001, ], opacity=0.3, color=(0, 0, 0.8))
		mlab.pipeline.iso_surface(src, contours=[-0.001, ], opacity=0.3, color=(0.8, 0, 0))
		mlab.show()
else:
	## Use matplotlib to show cuts of the molecular orbitals
	import matplotlib.pyplot as plt
	import numpy as np

	## Select cuts
	xd = mo_list[0][grid.N_[0]/2-1,:,:]
	yd = mo_list[0][:,grid.N_[1]/2-1,:]
	zd = mo_list[0][:,:,grid.N_[2]/2-1]

	## Plot cuts
	f, (pic1, pic2, pic3) = \
	    plt.subplots(3,1,sharex=True,sharey=True,figsize=(6,14))
	pic1.contour(z,y,xd,50,linewidths=0.5,colors='k')
	pic1.contourf(\
	    z,y,xd,50,cmap=plt.cm.rainbow,vmax=abs(xd).max(),vmin=-abs(xd).max())
	pic1.set_xlabel('z')
	pic1.set_ylabel('y')

	pic2.contour(z,x,yd,50,linewidths=0.5,colors='k')
	pic2.contourf(\
	    z,x,yd,50,cmap=plt.cm.rainbow,vmax=abs(yd).max(),vmin=-abs(yd).max())
	pic2.set_xlabel('z')
	pic2.set_ylabel('x')

	pic3.contour(y,x,zd,50,linewidths=0.5,colors='k')
	pic3.contourf(\
	    y,x,zd,50,cmap=plt.cm.rainbow,vmax=abs(zd).max(),vmin=-abs(zd).max())
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
