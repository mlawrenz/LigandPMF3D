LigandPMF3D
===========

** For now only used with MSMBuilder 2.7**

Python program for creating 3d PMFs and estimating binding free energies from protein-ligand binding simulations.

Program requires as input (see options with -h flag) the directiory with MSM output (Populations.dat, Mapping.dat), the location of the Gens.lh5 file, and the ligand index file.

Program output is a OpenDX file pmf3d_d1.dx (3D PMF at gridspace 1 Angstrom) centered on the protein as aligned in the Gens file.

