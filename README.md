# json-chem-tools

## Summary of files / Résumé des fichiers

* physics.py: Mainly unit conversion factors / Principalement des facteurs de conversion
* wavefunction.py: Function definitions for json_to_cube.py (pending review of ORBKIT) / Définitions de fonctions pour json_to_cube.py (attendant revue de ORBKIT)
* json_to_fchk.py: Script for converting scanlog .json file to Gaussian .fchk file / Script pour convertir un fichier .json de scanlog en .fchk de Gaussian
* json_to_cube.py: Script for discretising a wavefunction described in a scanlog .json file / Script pour discrétiser une fonction d'onde décrite dans un fichier .json de scanlog
* json_to_cube_orb.py: Same as json_to_cube.py, but based on ORBKIT / Identique à json_to_cube.py, mais basé sur ORBKIT
These two scripts, despite not having it explicit in their names, also process the data via ORBKIT:
* cclib_to_cube.py: Variant of json_to_cube.py accepting input via cclib / Variante de json_to_cube.py acceptant les données via cclib
* cclib_to_mayavi.py: Script for generating mayavi visualisation output from a chemistry file / Script pour générer une visualisation mayavi à partir d'un fichier de chimie
