#! /usr/bin/env python2.7
# -*- encoding: utf-8 -*-



## Usage
from sys import argv, stderr
if len(argv) < 3:
	stderr.write("Usage: {} $job[=$value] $npts < ${{input}}.json > ${{output}}.cube\n".format(argv[0]))
	exit(1)

from sys import stdin, stdout
from json import load
from math import pi, ceil
from physics import A_to_a0
from wavefunction import gs

## For Python 2 & 3 interoperability (`print` is a statement in Python 2)
write = stdout.write



## Get JSON data
struct = load(stdin)

molecule = struct["molecule"]
n_atoms = molecule["nb_atoms"]



## Size: gaussian cube protocol
temp = struct["results"]["geometry"]["elements_3D_coords_converged"]
x_max, y_max, z_max = max(temp[::3]), max(temp[1::3]), max(temp[2::3])
x_min, y_min, z_min = min(temp[::3]), min(temp[1::3]), min(temp[2::3])
t_x, t_y, t_z = x_max - x_min, y_max - y_min, z_max - z_min

## (p|s).
## Try to treat it as a number
try:
	p_npts = int(argv[2])
	## Two regimes: positive or zero, and negative
	## Positive or zero: the number of points is given and the spacing is deduced
	if 0 <= p_npts:
		if p_npts == 0:
			p_npts = 80
		p1, p2, p3 = (p_npts,)*3
		s_x, s_y, s_z = t_x/p1, t_y/p2, t_z/p3
		s1, s2, s3 = sorted([s_x, s_y, s_z])
	## Negative: the spacing is given and the number of points is deduced
	else:
		## -1 is not implemented
		s1, s2, s3 = (2.0**(2+p_npts)/3.0,)*3   if -5 < p_npts < -1 else \
		             (-p_npts*1e-3/A_to_a0,)*3  if p_npts <= -5 else None
		p1, p2, p3 = sorted([int(t_x/s1 + 1), int(t_y/s2 + 1), int(t_z/s3 + 1)])
	del p_npts

## Didn't work - it's a keyword
except ValueError:
	s1, s2, s3 = (1.0/3.0,)*3  if argv[2] == "Coarse" else \
	             (1.0/6.0,)*3  if argv[2] == "Medium" else \
	             (1.0/12.0,)*3 if argv[2] == "Fine" else None
	p1, p2, p3 = sorted([t_x/s1, t_y/s2, t_z/s3])



## Get job type and parameter
job, value = argv[1].split("=")

if job == "MO":
	## MO is always a list as this allows the same formatting to be used
	## for both single and multiple orbitals
	try:
		MO = [int(value)]
		n_val = 1
	except ValueError:
		if value == "All": pass
		elif value == "Valence": pass
		elif value == "Virtuals": pass
		elif value in {"OccA","OccB"}: pass
		elif value in {"AMO", "BMO"}: pass
		elif value in {"HOMO","LUMO"}:
			#def discrete(x, y, z):
			#	return sum(gs(l, m, a, x, y, z) for l in ... for m in ... for a in ...)
			MO = [(sum(molecule["atoms_Z"]) - molecule["charge"])//2 + (1 if value == "LUMO" else 0)]
			n_val = 1

elif job == "FDensity":
	#def discrete(x, y, z):
	#	return 
	n_val = 1

elif job == "Potential":
	pass



## Write title line
write(" {} {}\n".format(molecule["formula"], argv[1]))

## Write headers
write(" MO coefficients\n" if job == "MO" else
      " Electrostatic potential from Total SCF Density\n" if job == "Potential" and value == "SCF" else
      " Electron density from Total SCF Density\n" if job == "FDensity" and value == "SCF" else None)

write(" {:> 4d}   {:>9.6f}   {:>9.6f}   {:>9.6f}    {:d}\n".format(-n_atoms if job == "MO" else n_atoms, x_min, y_min, z_min, n_val))

line = " {:> 4d}   {:>9.6f}   {:>9.6f}   {:>9.6f}\n".format
write(line(p1, s1, 0, 0))
write(line(p2, 0, s2, 0))
write(line(p3, 0, 0, s3))

line = "  {:> 3d} {:>11.6f} {:>11.6f} {:>11.6f} {:>11.6f}\n".format
for a, p in zip(molecule["atoms_Z"], molecule["starting_geometry"]):
	write(line(a, a, p[0], p[1], p[2]))



## Specific to MO orbitals (c.f. cubegen documentation)
if job == "MO":
	write(" {:> 4d}".format(n_val))
	for o in MO:
		write("  {:> 3d}".format(o))
	write("\n")



## Generate voxel data
l = p3*n_val
block = ((" {: .5E}"*6 + "\n")*(l//6) + (" {: .5E}"*(l%6) + "\n" if l%6 > 0 else "")).format
write(block(*([0]*l)))

#X, Y, Z = [x_min + s_x*n for n in range(p1)], [y_min + s_y*n for n in range(p2)], [z_min + s_z*n for n in range(p3)]

#for x in X:
#	for y in Y:
#		for z in Z:
#			write(block(discrete(x, y, z)))