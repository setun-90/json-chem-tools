#! /usr/bin/env python2.7
## -*- encoding: utf-8 -*-



from physics import A_to_a0
import numpy as np

try:
	from enthought.mayavi import mlab
except ImportError as E:
	from sys import stderr
	stderr.write("import enthought.mayavi failed -- trying mayavi\n")
	try:
		from mayavi import mlab
	except ImportError:
		from sys import exit
		stderr.write("import mayavi failed\n")
		exit(1)



## Initialize visualization details common to all jobs
with open("Atoms.csv", "r") as f:
	tab = [line.split() for line in f]

def scene_init(j_data):
	figure = mlab.figure(bgcolor=(1,1,1))
	geom = np.array(j_data["results"]["geometry"]["elements_3D_coords_converged"]).reshape((-1,3))/A_to_a0

	if len(geom) > 1:
		## Eliminate hydrogens
		#mod_geom = geom[np.array(map(int, [p[1] for p in qc.geo_info if p[2] != '1.0']))]
		mod_geom = geom
		## Calculate best fitting plane via PCA
		eival, eivec = np.linalg.eig(np.cov((mod_geom - np.mean(mod_geom, axis=0)).T))
		sort = eival.argsort()
		eival, eivec = eival[sort], eivec[:,sort]
		print eival
		print eivec
		normal = eivec[:,0]
		## Grab point from best fitting plane (NOT the view) to use as focal point
		point = np.mean(geom, axis=0)

		from math import sqrt, atan2
		## Calculate hypotenuse on x-y plane
		hxy = sqrt(normal[0]**2 + normal[1]**2)
		## Calculate viewing distance r
		r = sqrt(hxy**2 + normal[2]**2)
		## Calculate azimuth a and elevation e
		## Python and Numpy use radians, but MayaVi uses degrees
		a, e = np.rad2deg(atan2(normal[0], normal[1])), np.rad2deg(atan2(r, hxy))
		mlab.view(azimuth=a, elevation=e-70, figure=figure, focalpoint=point)
		print a, e, point

		## DEBUG: show normal and view vectors
		#mlab.quiver3d(point[0], point[1], point[2], normal[0], normal[1], normal[2])
		#mlab.quiver3d([0]*3, [0]*3, [0]*3, [1,0,0], [0,1,0], [0,0,1], color=(0,0,0))
		#mlab.quiver3d(*np.concatenate((np.zeros((3,2)), eivec[1:])), color=(0,0,1))

	conn = j_data["molecule"]["connectivity"]["atom_pairs"]
	atom_nums = j_data["molecule"]["atoms_Z"]

	## Draw atoms and bonds
	for i, atom in enumerate(atom_nums):
		p, color = geom[i], tuple(float(x)/255.0 for x in tab[atom][3:6])

		## Requires >=MayaVi-4.6.0
		mlab.points3d([p[0]], [p[1]], [p[2]], figure=figure, mode='sphere', color=color, resolution=15, scale_factor=0.3)

	for pair in conn:
		att1 = tab[atom_nums[pair[0]]]
		p1, p2 = geom[pair[0]], geom[pair[1]]
		color = tuple(float(x)/255.0 for x in att1[3:6])

		mlab.quiver3d([p1[0]],         [p1[1]],         [p1[2]],
		              [p2[0] - p1[0]], [p2[1] - p1[1]], [p2[2] - p1[2]],
		              figure=figure, mode='cylinder', color=color, resolution=15, scale_factor=0.5)

	return figure



## Visualize
def topo(j_data, file_name=None):
	figure = scene_init(j_data)
	## Show labels and numbers ( = indices + 1 )
	#for i, atom in enumerate(j_data["molecule"]["atoms_Z"]):
	#	p, label = j_data["molecule"]["starting_geometry"][i], tab[atom][1]
	#	print p
	#	mlab.text3d(p[0], p[1], p[2], label, color=(0,0,0), figure=figure)

	if file_name is not None:
		e = mlab.view()[1]
		mlab.savefig("{}-TOPO.png".format(file_name), figure=figure)
	#	for i in range(10):
	#		mlab.view(elevation=e-i)
	#		mlab.savefig("{}-TOPO-{}.png".format(file_name, i))
	return figure

def viz_MO(data, X, Y, Z, j_data, file_name=None):
	figure = scene_init(j_data)
	for i, series in enumerate(data):

		MO_data = mlab.pipeline.scalar_field(X, Y, Z, series, figure=figure)

		MOp = mlab.pipeline.iso_surface(MO_data, figure=figure, contours=[ 0.05 ], color=(0.4, 0, 0.235))
		MOn = mlab.pipeline.iso_surface(MO_data, figure=figure, contours=[-0.05 ], color=(0.95, 0.90, 0.93))

		try:
			mlab.savefig("./{}-MO-{}.png".format(file_name, i))
		except:
			pass

		MOp.remove()
		MOn.remove()

	return figure

def viz_EDD(data, X, Y, Z, j_data, file_name=None):
	figure = scene_init(j_data)
	for i, series in enumerate(data):
		D_data = mlab.pipeline.scalar_field(X, Y, Z, series, figure=figure)

		Dp = mlab.pipeline.iso_surface(D_data, figure=figure, contours=[ 0.00025 ], color=(0.0, 0.5, 0.5))
		Dn = mlab.pipeline.iso_surface(D_data, figure=figure, contours=[-0.00025 ], color=(0.95, 0.95, 0.95))
		#Dn.actor.property.representation = 'wireframe'
		#Dn.actor.property.line_width = 0.5

		if file_name is not None:
			mlab.savefig("./{}-EDD-{}.png".format(file_name, i))

		Dp.remove()
		Dn.remove()

	return figure
