#! /usr/bin/env python2.7
## -*- encoding: utf-8 -*-



## Usage
from sys import argv, stderr
if len(argv) < 2:
	stderr.write("Usage: {} $job < ${{input}}.json > ${{output}}.fchk\n".format(argv[0]))
	exit(1)

from sys import stdin, stdout
from json import load
from physics import eV_to_Eh, cm_1_to_Eh
from wavefunction import primitives

## For Python 2 & 3 interoperability (`print` is a statement in Python 2)
write = stdout.write



## The main function for printing all lines in the output
#@profile   # for kernprof
def line(k, t, e):
	## t is necessary because certain arrays in the .json are
	## printed and formatted both as integers and as floats in
	## the .fchk
	## (e.g. "Integer atomic weights" and "Real atomic weights")

	## Try to treat it as an iterable
	try:
		l = len(e)

		if t == int:
			m = 6
			write(u"{:40}   I   N={:> 12d}\n".format(k, l))
			primitive = u" {: 11.0f}"

		elif t == float:
			m = 5
			write(u"{:40}   R   N={:> 12d}\n".format(k, l))
			primitive = u" {: .8E}"

		lr = l - l%m                                                    ## reduced length
		body = (primitive*m + "\n").format                              ## formatting function for a line
		tail = (primitive*(l - lr) + "\n" if lr < l else "").format     ## formatting function for remainder

		## Formatting and output are done line by line
		## to avoid duplicating the size of the array in
		## the format string (esp. the MO coefficients)
		for i in range(0, lr, m):
			write(body(*e[i:i+m]))
		write(tail(*e[lr:l]))

	## Didn't work - it's a scalar
	except TypeError:
		write(u"{:40}   I     {:> 12d}\n".format(k, e) if t == int else \
		      u"{:40}   R     {: 22.15E}\n".format(k, e))



## Cache certain properties
struct = load(stdin)
molecule = struct["molecule"]
results = struct["results"]
details = struct["comp_details"]



## Tests to check
def pred_a():
	return results["excited_states"]["et_energies"] != "N/A"

def is_closed_shell():
	return details["general"]["is_closed_shell"] == "True"

## $.fchk:1 (Name)
write(molecule["formula"] + "\n")

## $.fchk:2 (Information)
type = "FOPT"
write("{:10}{:30}                              {:30}\n".format(type, details["general"]["functional"], details["general"]["basis_set_name"]))

## The rest
n_atoms = molecule["nb_atoms"]
line("Number of atoms", int, n_atoms)
line("Info1-9", int, [9, 5, 0, 0, 0, 110, 1, 2, 18, -301])
line("Charge", int, molecule["charge"])
line("Multiplicity", int, molecule["multiplicity"])

n_electrons = sum(z for z in molecule["atoms_Z"]) - molecule["charge"]
n_alpha_electrons = (n_electrons + (molecule["multiplicity"] - 1)//2)//2
n_beta_electrons = (n_electrons - (molecule["multiplicity"] - 1)//2)//2
line("Number of electrons", int, n_electrons)
line("Number of alpha electrons", int, n_alpha_electrons)
line("Number of beta electrons", int, n_beta_electrons)

line("Number of basis functions", int, details["general"]["basis_set_size"])
line("Number of independent functions", int, results["wavefunction"]["MO_number"])
line("Number of point charges in /Mol/", int, 0)
line("Number of translation vectors", int, 0)
line("Atomic numbers", int, molecule["atoms_Z"])
line("Nuclear charges", float, molecule["atoms_Z"])
line("Current cartesian coordinates", float, results["geometry"]["elements_3D_coords_converged"])
line("Force Field", int, 0)
line("Int Atom types", int, [0]*n_atoms)
line("MM charges", float, [0]*n_atoms)
line("Integer atomic weights", int, molecule["atoms_masses"])
line("Real atomic weights", float, molecule["atoms_masses"])
line("Atom fragment info", int, [0]*molecule["nb_atoms"])
line("Atom residue num", int, [0]*molecule["nb_atoms"])
line("Nuclear spins", int, molecule["nuclear_spins"])
line("Nuclear ZEff", float, molecule["atoms_Zeff"])
line("Nuclear ZNuc", float, molecule["atoms_Z"])
line("Nuclear QMom", float, molecule["nuclear_QMom"])
line("Nuclear GFac", float, molecule["nuclear_gfactors"])

line("MicOpt", int, [-1]*n_atoms)

## Extract contracted shells (both S=P and normal)
primitives_by_shell = primitives(details["general"]["basis_set"])

## Since the extraction is just an optimization (cclib does not coalesce hybridized orbitals),
## leave this line for discretisations not taking this into account (like ours!)
#primitives_by_shell = [(i, ssh) for i, atom in enumerate(details["general"]["basis_set"], 1) for ssh in atom]
shell_type = {
	"S"  : 0,
	"P"  : 1,
	"SP" :-1,
	"D"  :-2,
	"F"  :-3,
	"G"  :-4
}
shell_types = [shell_type[ssh[1][0]] for ssh in primitives_by_shell]
primitives_per_shell = [len(p[1][1]) for p in primitives_by_shell]
shell_to_atom = [atom[0] for atom in primitives_by_shell]

line("Number of contracted shells", int, len(primitives_by_shell))
line("Number of primitive shells", int, sum(len(prim[1][1]) for prim in primitives_by_shell))

line("Pure/Cartesian d shells", int, 0)
line("Pure/Cartesian f shells", int, 0)

line("Highest angular momentum", int, max(map(abs, shell_types)))
line("Largest degree of contraction", int, max(primitives_per_shell))
line("Shell types", int, shell_types)
line("Number of primitives per shell", int, primitives_per_shell)
line("Shell to atom map", int, shell_to_atom)
line("Primitive exponents", float, [C[0] for prim in primitives_by_shell for C in prim[1][1]])
line("Contraction coefficients", float, [C[1] for prim in primitives_by_shell for C in prim[1][1]])
line("P(S=P) Contraction coefficients", float, [C[2] for prim in primitives_by_shell for C in prim[1][1]])
line("Coordinates of each shell", float, [e for i in shell_to_atom for e in results["geometry"]["elements_3D_coords_converged"][3*i-3:3*i]])

## To be revised later
line("Num ILSW", int, 100)
## Massive constant array that varies in only three entries
line("ILSW", int, [0, 	    0,       0,       0,       2,       0]
                + [0, 	    0,       0,       0,    1009,      -1]
                + [0,       0,       0,       0,       2,       0]    # some jobs have 0 instead of 2
                + [0, 	    0,       0,       0,       1,       0]    # some jobs have 0 instead of 1
                + [1,       0,       0,       0,       0,       0]
                + [0,       0,   10000,       0,      -1,       0]
                + [0,       0,       0,       0,       0,       0]
                + [0,       0,       0,       1,       0,       0]
                + [0,       0,       1,       0,       0,       0]
                + [0,       0,       4,      41,       0,       0]
                + [0,       0,      13,       0,       0,       0]
                + [0,       0,       0, n_atoms,       0,       0]
                + [0,       0,       0,       0,       0,       0]*4
                + [0,       0,       0,       0])

line("Num RLSW", int, 41)
line("RLSW", float, [.75]*2 + [1]*2 + [.25]
                  + ([0]*2 + [1] + [0]*2)*6
                  + [0]*3 + [1]*2 + [0])

mx_bond = max(molecule["atoms_valence"])
if mx_bond == 0:
	mx_bond = 1

line("MxBond", int, mx_bond)
line("NBond", int, molecule["atoms_valence"])

## DO NOT CHANGE TO `[[0]*mx_bond]*n_atoms` (`*` creates references not copies)
i_bond = [[0   for i in range(mx_bond)] for i in range(n_atoms)]
r_bond = [[0.0 for i in range(mx_bond)] for i in range(n_atoms)]
ind = [0]*n_atoms
for p_atom, order in zip(molecule["connectivity"]["atom_pairs"], molecule["connectivity"]["bond_orders"]):
	i = p_atom[0]
	i_bond[i][ind[i]] = p_atom[1] + 1
	r_bond[i][ind[i]] = order
	ind[i] += 1

line("IBond", int, [e2 for e1 in i_bond for e2 in e1])
line("RBond", float, [e2 for e1 in r_bond for e2 in e1])
line("Virial Ratio", float, results["wavefunction"]["virial_ratio"][-1])

if results["wavefunction"]["total_molecular_energy"] != "N/A":
	SCF_energy = results["wavefunction"]["total_molecular_energy"]/eV_to_Eh
	line("SCF Energy", float, SCF_energy)
	if pred_a():
		CIS_energy = SCF_energy + results["excited_states"]["et_energies"][0]/cm_1_to_Eh
		line("CIS Energy", float, CIS_energy)
		line("Total Energy", float, CIS_energy)
	#	line("Post-SCF wavefunction norm", float, )          # should be scalar
	else:
		line("Total Energy", float, SCF_energy)

if results["geometry"]["geometric_values"][0][0] != "N/A":
	line("RMS Force", float, results["geometry"]["geometric_values"][-1][1])
line("RMS Density", float, results["wavefunction"]["SCF_values"][0])

line("External E-field", float, [0]*35)
line("IOpCl", int, 0 if is_closed_shell() else 1)

mo_coeffs = results["wavefunction"]["MO_coefs"]
mo_energies = results["wavefunction"]["MO_energies"]
## Check condition: 1 if n_alpha_electrons != n_beta_electrons and len(mo_energies) == 1 else 0 ???
line("IROHF", int, 1 if n_alpha_electrons != n_beta_electrons and len(mo_energies) == 1 else 0)
line("Alpha Orbital Energies", float, [E/eV_to_Eh for E in mo_energies[0]])
if not is_closed_shell():
	line("Beta Orbital Energies", float, [E/eV_to_Eh for E in mo_energies[1]])
line("Alpha MO coefficients", float, [C for arr in mo_coeffs for C in arr])
#if argv[1] in {"SP_moins", "SP_plus"}:
#	line("Beta MO coefficients", float, [C for arr in mo_coeffs for C in arr])
l = details["general"]["basis_set_size"]
line("Total SCF Density", float, [2*sum(mo_coeffs[i][a]*mo_coeffs[i][b] for i in range(n_electrons//2)) for a in range(l) for b in range(a,l)])
#line("Total Cl RHo(1) Density", float, )
#line("Total Cl Density", float, )

if results["wavefunction"]["Mulliken_partial_charges"] != "N/A":
	line("Mulliken Charges", float, results["wavefunction"]["Mulliken_partial_charges"])

line("Optimization MaxStp", int, 100)
line("Optimization Job offset", int, 0)
line("Optimization Num results per geometry", int, 2)
line("Optimization Num geometry variables", int, n_atoms*3)
#line("Opt point       1 Results for each geome", float, )
#line("Opt point       1", float,)
#line("Optimization Number of geometries", int, )
line("ONIOM Charges", int, [0]*16)
line("ONIOM Multiplicities", int, [0]*16)
line("Atom Layers", int, [1]*n_atoms)
line("Atom Modifiers", int, [0]*n_atoms)
line("Force Field", int, 0)
line("Int Atom Modified Types", int, [0]*n_atoms)
line("Link Atoms", int, [0]*n_atoms)
line("Atom Modified MM Charges", float, [0]*n_atoms)
line("Link Distances", float, [0]*n_atoms)
#line("Cartesian Gradient", float, )

#if argv[1] in {"OPT", "SP_plus", "SP_moins", "OPT_ES", "OPT_ET"}:
#	line("Dipole moment", float, )
#	line("Quadrupole moment", float, )
#	line("QEq coupling tensors", float, )

#if is_closed_shell():
#	line("Anisotropic Hyperfine tensors", float, )
#	line("Isotropic Hyperfine splittings", float, )

if pred_a():
#	line("ETran NETS", int, details["excited_states"]["nb_et_states"])
#	line("ETran NETV", int, )
#	line("ETran scalars", int, )
#	line("ETran spin", int, )
	line("ETran sym", int, results["excited_states"]["et_sym"])
#	line("ETran state values", float, )
