#! /usr/bin/env python2.7
## -*- encoding: utf-8 -*-



## Usage
from sys import argv, stderr
if len(argv) < 5:
	stderr.write("Usage: {} $job[=$value] $npts ${{input}}.log ${{output}}.cube\n".format(argv[0]))
	exit(1)

from sys import stdin, stdout
from math import pi, ceil
from physics import A_to_a0
from orbkit import options, grid, read, output, core

## For Python 2 & 3 interoperability (`print` is a statement in Python 2)
write = stdout.write



## Set input/output
from cclib.parser.gaussianparser import Gaussian
data = Gaussian(argv[3]).parse()
qc = read.convert_cclib(data)



## Size: gaussian cube protocol
over_s = 5   ## to be tuned

## (p|s).
## Try to treat it as a number
try:
	p_npts = int(argv[2])
	## Two regimes: positive or zero, and negative
	## Positive or zero: the number of points is given and the spacing is deduced
	if 0 <= p_npts:
		if p_npts == 0:
			p_npts = 80
		## ORBKIT
		grid.N_ = [(p_npts,)*3]

	## Negative: the spacing is given and the number of points is deduced
	else:
		## ORBKIT
		## -1 is not implemented
		grid.adjust_to_geo(qc, over_s, 2.0**(2+p_npts)/3.0   if -5 < p_npts < -1 else \
		                               -p_npts*1e-3/A_to_a0  if p_npts <= -5 else None)

	del p_npts

## Didn't work - it's a keyword
except ValueError:
	## ORBKIT
	options.adjust_grid = [over_s, 1.0/3.0  if argv[2] == "Coarse" else \
	                               1.0/6.0  if argv[2] == "Medium" else \
	                               1.0/12.0 if argv[2] == "Fine" else None]



## Get job type and parameter
job, value = argv[1].split("=")

if job == "MO":
	## Get list of orbitals
	try:
		options.calc_mo = map(int, value.split(","))
	except ValueError:
		if value == "All":
			options.calc_mo = ["all_mo"]
		elif value in {"HOMO","LUMO"}:
			options.calc_mo = [value.lower()]

elif job == "FDensity":
	pass

elif job == "Potential":
	pass

out = core.rho_compute(qc, numproc=4)

output.cube_creator(out, argv[4], qc.geo_info, qc.geo_spec)
