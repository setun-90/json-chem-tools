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
from numpy import ndarray
from orbkit import options, grid, read, output, core

## For Python 2 & 3 interoperability (`print` is a statement in Python 2)
write = stdout.write



## Set input/output
from cclib.parser import ccopen
data = ccopen(argv[3]).parse()
qc = read.convert_cclib(ccopen(argv[3]).parse(), all_mo=True)



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
		grid.N_ = [(p_npts,)*3]
		#grid.

	## Negative: the spacing is given and the number of points is deduced
	else:
		## -1 is not implemented
		grid.adjust_to_geo(qc, over_s, 2.0**(2+p_npts)/3.0   if -5 < p_npts < -1 else \
		                               -p_npts*1e-3/A_to_a0  if p_npts <= -5 else None)

	del p_npts

## Didn't work - it's a keyword
except ValueError:
	grid.adjust_to_geo(qc, over_s, 1.0/3.0  if argv[2] == "Coarse" else \
	                               1.0/6.0  if argv[2] == "Medium" else \
	                               1.0/12.0 if argv[2] == "Fine" else None)



## Get job type and parameter
job, value = argv[1].split("=")

if job == "MO":
	## Get list of orbitals
	try:
		qc.mo_spec = read.mo_select(qc.mo_spec, map(int, value.split(",")))["mo_spec"]
	except ValueError:
		if value == "All":
			qc.mo_spec = read.mo_select(qc.mo_spec, ["all_mo"])["mo_spec"]
		elif value in {"HOMO","LUMO"}:
			qc.mo_spec = read.mo_select(qc.mo_spec, [value.lower()])["mo_spec"]
	def func(data):
		return core.rho_compute(qc, calc_mo=True, numproc=4)

elif job == "FDensity":
	def func(data):
		return core.rho_compute(qc, numproc=4)

elif job == "Potential":
	pass

out = func(qc)

print out[0]
output.main_output(out[0,:], qc.geo_info, qc.geo_spec, outputname=argv[4], otype='cb')
