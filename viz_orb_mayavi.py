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

## Elevation angle
angle = 20

def _scene_init(j_data):
	figure = mlab.figure(bgcolor=(1,1,1))
	figure.scene.disable_render= True
	geom = np.array(j_data["results"]["geometry"]["elements_3D_coords_converged"]).reshape((-1,3))/A_to_a0

	if len(geom) > 1:
		## Eliminate hydrogens
		#mod_geom = geom[np.array(map(int, [p[1] for p in qc.geo_info if p[2] != '1.0']))]
		mod_geom = geom
		## Calculate best fitting plane via PCA
		eival, eivec = np.linalg.eig(np.cov((mod_geom - np.mean(mod_geom, axis=0)).T))
		sort = eival.argsort()
		eival, eivec = eival[sort], eivec[:,sort]
		normal = eivec[:,0]
		## Grab point from best fitting plane (NOT the view) to use as focal point
		#point = np.mean(geom, axis=0)

		from math import sqrt, acos, atan2
		## Calculate viewing distance r
		r = sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2)
		## Calculate azimuth a and elevation e
		## Python and Numpy use radians, but MayaVi uses degrees
		a, e = np.rad2deg(atan2(normal[1], normal[0])), np.rad2deg(acos(normal[2]/r))
		mlab.view(azimuth=a, elevation=e+angle, figure=figure)

		## DEBUG: show normal and view vectors
		#print a, e
		#print mlab.view()
		#mlab.quiver3d(point[0], point[1], point[2], normal[0], normal[1], normal[2])
		#mlab.quiver3d([0]*3, [0]*3, [0]*3, [1,0,0], [0,1,0], [0,0,1], color=(0,0,0))
		#mlab.quiver3d(*np.concatenate((np.zeros((3,2)), eivec[1:])), color=(0,0,1))

	else:
		normal = np.array([0,0,1])
		mlab.view(azimuth=0, elevation=0, figure=figure)

	conn = j_data["molecule"]["connectivity"]["atom_pairs"]
	atom_nums = j_data["molecule"]["atoms_Z"]

	## Draw atoms and bonds
	for i, atom in enumerate(atom_nums):
		p, color = geom[i], tuple(float(x)/255.0 for x in tab[atom][3:6])

		## Requires >=MayaVi-4.6.0
		mlab.points3d([p[0]], [p[1]], [p[2]],
		              figure=figure, mode='sphere', color=color, resolution=15, scale_factor=0.5)

	for pair in conn:
		att1 = tab[atom_nums[pair[0]]]
		p1, p2 = geom[pair[0]], geom[pair[1]]
		color = tuple(float(x)/255.0 for x in att1[3:6])

		mlab.quiver3d([p1[0]],         [p1[1]],         [p1[2]],
		              [p2[0] - p1[0]], [p2[1] - p1[1]], [p2[2] - p1[2]],
		              figure=figure, mode='cylinder', color=color, resolution=15, scale_factor=0.5)

	#mlab.axes(figure=figure)

	return figure, normal



## Visualize
def topo(j_data, file_name=None, size=(600,600)):
	"""Creates the topological view of the molecule.
	** Parameters **
	  j_data : dict
	Data on the molecule, as deserialized from the scanlog format.
	  file_name : str, optional
	Base name of the file in which to save the image.
	  size : tuple(int, int), optional
	The size of the image to save.
	** Returns **
	  figure
	The MayaVi scene containing the visualization.
	"""

	figure, normal = _scene_init(j_data)
	geom = np.array(j_data["results"]["geometry"]["elements_3D_coords_converged"]).reshape((-1,3))/A_to_a0

	## Show labels and numbers ( = indices + 1 )
	for i, atom in enumerate(j_data["molecule"]["atoms_Z"]):
		P, label = geom[i], tab[atom][1]
		mlab.text3d(P[0] - normal[0], P[1] - normal[1], P[2] - normal[2], label + str(i + 1), color=(0,0,0), scale=0.5, figure=figure)

	if file_name is not None:
		mlab.savefig("{}-TOPO.png".format(file_name), figure=figure, size=size)

	return figure

def viz_MO(data, X, Y, Z, j_data, file_name=None, labels=None, size=(600,600)):
	"""Visualizes the molecular orbitals of the molecule.
	** Parameters **
	  data : list(numpy.ndarray)
	List of series of voxels containing the scalar values of the molecular orbitals to plot.
	  X, Y, Z
	Meshgrids as generated by numpy.mgrid, for positioning the voxels.
	  j_data : dict
	Data on the molecule, as deserialized from the scanlog format.
	  file_name : str, optional
	Base name of the files in which to save the images.
	  labels : list(str), optional
	Labels to append to `file_name`. If None, the index of the series is appended.
	  size : tuple(int, int), optional
	The size of the image to save.
	** Returns **
	  figure
	The MayaVi scene containing the visualization.
	"""

	figure = _scene_init(j_data)[0]
	for i, series in enumerate(data):

		MO_data = mlab.pipeline.scalar_field(X, Y, Z, series, figure=figure)

		MOp = mlab.pipeline.iso_surface(MO_data, figure=figure, contours=[ 0.05 ], color=(0.4, 0, 0.235))
		MOn = mlab.pipeline.iso_surface(MO_data, figure=figure, contours=[-0.05 ], color=(0.95, 0.90, 0.93))

		if file_name is not None:
			mlab.savefig("./{}-MO-{}.png".format(file_name, labels[i] if labels is not None else i), figure=figure, size=size)

		MOp.remove()
		MOn.remove()

	return figure

def viz_EDD(data, X, Y, Z, j_data, file_name=None, labels=None, size=(600,600)):
	"""Visualizes the electron density differences for the transitions of the molecule.
	** Parameters **
	  data : list(numpy.ndarray)
	Voxels containing the scalar values of the electron density differences to plot.
	  X, Y, Z
	Meshgrids as generated by numpy.mgrid, for positioning the voxels.
	  j_data : dict
	Data on the molecule, as deserialized from the scanlog format.
	  file_name : str, optional
	Base name of the files in which to save the images.
	  labels : list(str), optional
	Labels to append to `file_name`. If None, the index of the series is appended.
	  size : tuple(int, int), optional
	The size of the image to save.
	** Returns **
	  figure
	The MayaVi scene containing the visualization.
	"""

	figure = _scene_init(j_data)[0]
	for i, series in enumerate(data):
		D_data = mlab.pipeline.scalar_field(X, Y, Z, series, figure=figure)

		Dp = mlab.pipeline.iso_surface(D_data, figure=figure, contours=[ 0.00025 ], color=(0.0, 0.5, 0.5))
		Dn = mlab.pipeline.iso_surface(D_data, figure=figure, contours=[-0.00025 ], color=(0.95, 0.95, 0.95))
		#Dn.actor.property.representation = 'wireframe'
		#Dn.actor.property.line_width = 0.5

		if file_name is not None:
			mlab.savefig("./{}-EDD-{}.png".format(file_name, labels[i] if labels is not None else i), figure=figure, size=size)

		Dp.remove()
		Dn.remove()

	return figure

def viz_BARY(data, j_data, file_name=None, labels=None, size=(600,600)):
	"""Visualizes the barycenters of the electron density difference (for visualizing dipole moments).
	** Parameters **
	  data : tuple(numpy.ndarray((3,N)), numpy.ndarray((3,N)))
	Pair of column-major matrices containing the coordinates of the positive and negative barycenters, in that order.
	  j_data : dict
	Data on the molecule, as deserialized from the scanlog format.
	  file_name : str, optional
	Base name of the files in which to save the images.
	  labels : list(str), optional
	Labels to append to `file_name`. If None, the index of the series is appended.
	  size : tuple(int, int), optional
	The size of the image to save.
	** Returns **
	  figure
	The MayaVi scene containing the visualization.
	"""

	figure = _scene_init(j_data)[0]

	for i, D in enumerate(data):
		Pp = mlab.points3d(D[0][0], D[0][1], D[0][2], figure=figure, color=(0.0, 0.5, 0.5))
		Pm = mlab.points3d(D[1][0], D[1][1], D[1][2], figure=figure, color=(0.95, 0.95, 0.95))

		if file_name is not None:
			mlab.savefig("./{}-BARY-{}.png".format(file_name, labels[i] if labels is not None else i), figure=figure, size=size)

		Pp.remove()
		Pm.remove()

	return figure
