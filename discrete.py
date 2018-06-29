#! /usr/bin/env python2.7
## -*- encoding: utf-8 -*-



from sys import argv
from physics import A_to_a0
from calc_orb import *
from viz_orb_mayavi import *
## Patched version of orbkit.read
import numpy as np



## Usage
usage = ( "Usage: {0} $action MO[=$value]                  $npts ${{input}}.json                      ${{output}}\n"
        + "     | {0} $action Potential                    $npts ${{input}}.json                      ${{output}}\n"
        + "     | {0} $action (TD|EDD|BARY|Tozer)[=$value] $npts ${{OPT-input}}.json ${{TD-input}}.json ${{output}}\n"
	+ "     | {0} $action topo                         $npts ${{input}}.json                      ${{output}}\n").format(argv[0])

if len(argv) < 5:
	from sys import stderr, exit
	stderr.write(usage)
	exit(0)



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



## Get grid parameters and initialize grid

## Oversizing in Bohr radii
## ORBKIT uses 5 by default, tune this as required
over_s = 7

## Spacing/Number of points
par = int(argv[3])



## Get job type and parameters
val = argv[2].split("=")
job = val[0]



## Calculate
if job == "topo":
	try:
		file_name=argv[5]
	except IndexError:
		file_name=None

	topo(data, file_name)

elif job == "MO":
	## Get list of orbitals
	MO_list = val[1].split(",")

	out, X, Y, Z = MO(data, MO_list, grid_par=par)
	for series in out:
		## The length product works because all voxels of the ORBKIT grid have the same dimensions
		print np.sum(np.square(series))*(X[1,0,0] - X[0,0,0])*(Y[0,1,0] - Y[0,0,0])*(Z[0,0,1] - Z[0,0,0])

	try:
		file_name=argv[5]
	except IndexError:
		file_name=None

	viz_MO(out, X, Y, Z, data, file_name=file_name)

elif job in {"TD", "EDD", "BARY", "Tozer"}:
	if ".json" not in argv[5]:
		from sys import stderr, exit
		stderr.write(usage)
		exit(1)

	with open(argv[5], "r") as f:
		TD_data = load(f)

	try:
		T_list = map(int, val[1].split(","))
		transitions = [TD_data["results"]["excited_states"]["et_transitions"][i] for i in T_list]
	except IndexError:
		transitions = TD_data["results"]["excited_states"]["et_transitions"]


	try:
		file_name=argv[6]
	except IndexError:
		file_name=None

	out, X, Y, Z = TD(data, transitions, grid_par=par)

	if job == "TD":
		for l in [(e[1], e[2][:3]) for e in out]:
			print l

	elif job == "EDD":
		viz_EDD([e[0] for e in out], X, Y, Z, data, file_name=file_name)

	elif job == "BARY":
		B = [e[2][3:] for e in out]
		for P in B:
			print P
		viz_BARY(B, data, file_name=file_name)

	elif job == "Tozer":
		print [e[1] for e in out]

elif job == "Potential":
	out_r, out_V, X, Y, Z = Potential(data, grid_par=par)

	viz_Potential(out_V, X, Y, Z, data)
