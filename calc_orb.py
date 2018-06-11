#! /usr/bin/env python2.7
## -*- encoding: utf-8 -*-



from physics import A_to_a0
from orbkit import grid, core
## Patched version of orbkit.read
import read
import numpy as np



## Get grid parameters and initialize grid
def init_ORB_grid(data, grid_par=-6, over_s=7):

	## Spacing/Number of points
	if grid_par > 0:
		grid.N_ = [grid_par]*3
	elif grid_par == 0:
		grid.N_ = [80]*3
	else:
		grid.delta_ = [1.0/(-grid_par)]*3

	grid.max_ = np.amax(data.geo_spec.T, axis=1) + over_s
	grid.min_ = np.amin(data.geo_spec.T, axis=1) - over_s

	grid.init()

	return np.mgrid[grid.min_[0]:grid.max_[0]:1j*len(grid.x),   # X
                        grid.min_[1]:grid.max_[1]:1j*len(grid.y),   # Y
                        grid.min_[2]:grid.max_[2]:1j*len(grid.z)]   # Z



## Calculations
def MO(j_data, MO_list, grid_par=-6):
	qc = read.convert_json(j_data, all_mo=True)
	X, Y, Z = init_ORB_grid(qc, grid_par=grid_par)

	## Get list of orbitals
	qc.mo_spec = read.mo_select(qc.mo_spec, MO_list)["mo_spec"]

	## Calculate
	out = core.rho_compute(qc, calc_mo=True, numproc=4)

	return out, X, Y, Z



def EDD(j_data, T_data, T_list, grid_par=-6):
	transitions = T_data["results"]["excited_states"]["et_transitions"]

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
	X, Y, Z = init_ORB_grid(qc, grid_par=grid_par)
	MOs = core.rho_compute(qc, calc_mo=True, numproc=4)

	## 2. Combine MOs according to info in `et_transitions`
	out = []
	for i, T in enumerate(transitions):
		series = np.zeros(MOs[0].shape)
		for j, ST in enumerate(T):
			## Dp_i = S_j(C_ij**2*(MO2_ij**2 - MO1_ij**2))
			print "Calculating transition {}.{}".format(i, j)
			series += ST[2]**2*(np.square(MOs[tab[ST[1][0]]]) - np.square(MOs[tab[ST[0][0]]]))
		out += [series]

	return out, X, Y, Z
