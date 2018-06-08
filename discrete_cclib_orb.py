#! /usr/bin/env python2.7
## -*- encoding: utf-8 -*-



from sys import argv
from physics import A_to_a0
from orbkit import grid, core, output, extras
## Patched version of orbkit.read
import read
import numpy as np



## Usage
usage = ( "Usage: {0} $action MO[=$value]  $npts ${{input}}.(log|json) ${{output}}\n"
        + "     | {0} $action EDD[=$value] $npts ${{OPT-input}}.json   ${{TD-input}}.json ${{output}}\n"
	+ "     | {0} $action topo         $npts ${{input}}.json ${{output}}\n").format(argv[0])



## Input
from json import load
if ".json" not in argv[4]:
	from sys import stderr, exit
	stderr.write(usage)
	exit(1)

	## Code of uncertain status due to instability in specs
	## Kept for later
	#from cclib.parser import ccopen
	#cc_data = ccopen(argv[4]).parse()
	##cc_data.mocoeffs[0] = [[0.0 if np.log(np.abs(x)) < -9 else x for x in S] for S in cc_data.mocoeffs[0]]
	#qc = read.convert_cclib(cc_data, all_mo=True)

with open(argv[4], "r") as f:
	data = load(f)
qc = read.convert_json(data, all_mo=True)



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

grid.max_ = np.amax(qc.geo_spec.T, axis=1) + over_s
grid.min_ = np.amin(qc.geo_spec.T, axis=1) - over_s

grid.init()



## Initialize visualization details common to all jobs
try:
	from enthought.mayavi import mlab
except ImportError as E:
	from sys import stderr
	stderr.write("import enthought.mayavi failed -- trying mayavi")
	try:
		from mayavi import mlab
	except ImportError:
		from sys import exit
		stderr.write("import mayavi failed")
		exit(1)

mlab.figure(bgcolor=(1,1,1))

if len(qc.geo_spec) > 1:
	## Eliminate hydrogens
	mod_geo_spec = qc.geo_spec[np.array(map(int, [p[1] for p in qc.geo_info if p[2] != '1.0']))]
	## Calculate best fitting plane via PCA
	eival, eivec = np.linalg.eig(np.cov((qc.geo_spec - np.mean(qc.geo_spec, axis=0)).T))
	sort = eival.argsort()[::-1]
	eival, eivec = eival[sort], eivec[:,sort]
	normal = eivec[:,-1]

	from math import sqrt, atan2
	## Calculate projection of normal on best fitting plane
	r_p = sqrt(normal[0]**2 + normal[1]**2)
	## Calculate viewing distance r
	r = sqrt(r_p**2 + normal[2]**2)
	## Calculate azimuth a and elevation e
	a, e = atan2(normal[0], normal[1]), atan2(r, r_p)
	mlab.view(azimuth=a, elevation=e-20, distance=min(grid.max_))

conn = data["molecule"]["connectivity"]["atom_pairs"]
atom_nums = data["molecule"]["atoms_Z"]
geom = data["molecule"]["starting_geometry"]
with open("Atoms.csv", "r") as f:
	tab = [line.split() for line in f]

## Draw atoms and bonds
for i, atom in enumerate(atom_nums):
	p, color = geom[i], tuple(float(x)/255.0 for x in tab[atom][3:6])

	mlab.points3d([p[0]], [p[1]], [p[2]], mode='sphere', color=color, resolution=15, scale_factor=0.3)

for pair in conn:
	att1 = tab[atom_nums[pair[0]]]
	p1, p2 = geom[pair[0]], geom[pair[1]]
	color = tuple(float(x)/255.0 for x in att1[3:6])

	mlab.quiver3d([p1[0]],         [p1[1]],         [p1[2]],
	              [p2[0] - p1[0]], [p2[1] - p1[1]], [p2[2] - p1[2]],
	              mode='cylinder', color=color, resolution=15, scale_factor=0.5)

X, Y, Z = np.mgrid[grid.min_[0]:grid.max_[0]:1j*len(grid.x),
                   grid.min_[1]:grid.max_[1]:1j*len(grid.y),
                   grid.min_[2]:grid.max_[2]:1j*len(grid.z)]



## Get job type and parameters
val = argv[2].split("=")
job = val[0]



## Calculate
if job == "topo":
	## Show labels and numbers ( = indices + 1 )
	for i, atom in enumerate(atom_nums):
		p, label = geom[i], tab[atom][1]
		mlab.text3d([p[0]], [p[1]], [p[2]], label, color=(0,0,0))

	mlab.savefig("{}-TOPO.png".format(argv[5]))

elif job == "MO":
	## Get list of orbitals
	qc.mo_spec = read.mo_select(qc.mo_spec, val[1].split(","))["mo_spec"]

	out = core.rho_compute(qc, calc_mo=True, numproc=4)
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
		for i, series in enumerate(out):

			MO_data = mlab.pipeline.scalar_field(X, Y, Z, series)

			MOp = mlab.pipeline.iso_surface(MO_data, contours=[ 0.05 ], color=(0.4, 0, 0.235))
			MOn = mlab.pipeline.iso_surface(MO_data, contours=[-0.05 ], color=(0.95, 0.90, 0.93))

			mlab.savefig("./{}-MO-{}.png".format(argv[5], i))
			MOp.remove()
			MOn.remove()

elif job == "EDD":
	try:
		T_list = val[1].split(",")
	except IndexError:
		T_list = ["*"]

	if ".json" not in argv[5]:
		from sys import stderr, exit
		stderr.write(usage)
		exit(1)

	with open(argv[5], "r") as f:
		TD_data = load(f)

	transitions = TD_data["results"]["excited_states"]["et_transitions"]

	## To save time, the calculation is done in two phases:

	## 1. Get all MOs involved in transitions and calculate them once
	MO_list = [STT[0] for T in transitions for ST in T for STT in ST[:2]]
	MO_set = set(MO_list)
	## The dictionary is needed because the MO numbers do not correspond
	## to their indices in the output of rho_compute, and we can only access
	## them through indices (assuming ORBKIT conserves the order of the MO numbers
	## in qc.mo_spec)
	tab = dict([(MO, i) for i, MO in enumerate(MO_set)])
	qc.mo_spec = read.mo_select(qc.mo_spec, [str(MO + 1) for MO in MO_set])["mo_spec"]
	MOs = core.rho_compute(qc, calc_mo=True, numproc=4)

	## 2. Combine MOs according to info in `et_transitions`
	for i, T in enumerate(transitions):
		series = np.zeros(MOs[0].shape)
		for j, ST in enumerate(T):
			## Dp_i = S_j(C_ij**2*(MO2_ij**2 - MO1_ij**2))
			print "Calculating transition {}.{}".format(i, j)
			series += ST[2]**2*(np.square(MOs[tab[ST[1][0]]]) - np.square(MOs[tab[ST[0][0]]]))

		D_data = mlab.pipeline.scalar_field(X, Y, Z, series)

		Dp = mlab.pipeline.iso_surface(D_data, contours=[ 0.00025 ], color=(0.0, 0.5, 0.5))
		Dn = mlab.pipeline.iso_surface(D_data, contours=[-0.00025 ], color=(0.95, 0.95, 0.95))
		#Dn.actor.property.representation = 'wireframe'
		#Dn.actor.property.line_width = 0.5

		mlab.savefig("./{}-EDD-{}.png".format(argv[6], i))
		Dp.remove()
		Dn.remove()
