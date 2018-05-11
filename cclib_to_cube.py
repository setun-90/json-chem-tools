#! /usr/bin/env python2.7
## -*- encoding: utf-8 -*-



## Usage
from sys import argv
if len(argv) < 5:
	from sys import stderr, exit
	stderr.write("Usage: {} $job[=$value] $npts ${{input}}.log ${{output}}\n".format(argv[0]))
	exit(1)

from physics import A_to_a0
from orbkit import grid, core, read, output, extras



## Input
from cclib.parser import ccopen
qc = read.convert_cclib(ccopen(argv[3]).parse(), all_mo=True)



## Get grid parameters and initialize grid
over_s = 5   # to be tuned

## Try to treat size parameter as a number
try:
	p_npts = int(argv[2])
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

grid.max_ = map(lambda a: max(a) + over_s, qc.geo_spec)
grid.min_ = map(lambda a: min(a) - over_s, qc.geo_spec)
grid.init()



## Get job type and parameters
job, value = argv[1].split("=")

if job == "MO":
	## Get list of orbitals
	try:
		mos = map(int, value.split(","))
	except ValueError:
		if value == "All":
			mos = ["all_mo"]
		elif value in {"HOMO","LUMO"}:
			mos = [value.lower()]
	#qc.mo_spec = read.mo_select(qc.mo_spec, mos)["mo_spec"]
	def func(data):
		#return core.rho_compute(data, calc_mo=True, numproc=4)
		return extras.calc_mo(data, mos, numproc=4)[0]

elif job == "FDensity":
	def func(data):
		return core.rho_compute(data, numproc=4)

elif job == "Potential":
	pass



## Main calculation
out = func(qc)



## Output
output.main_output(out[0,:], qc.geo_info, qc.geo_spec, outputname=argv[4], otype="cb")
