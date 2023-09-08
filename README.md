Material for paper: 
A continuum model for the elongation and orientation of Von
Willebrand Factor with applications in arterial flow

Link to archive article: 

Author of files in this repository: E. Yeo

Contents:

Python files (requires fenics2019.2.0.dev0 which is imported through dolfin):

    main.py

GMSH files (requires gmsh-3.0.6):


MATLAB Files (requires MATLAB lisence) 3) param.m 4) plots.m 5) Data

Raw Data exported python:



    main file runs all simulations for data presented in paper. Each figure is produced by a different 'study' appropriately labeled. Data for each figure is plotted in the results section. Data was then exported to MATLAB to provide uniform formatting and visualisation.

    Numerical verification file to suport main file. Model unchanged, mesh convergence is tested through simulations with increasingly fine meshes. Results sensitivity to the approximate parameter 'delta' is tested.

    Contains parameters used to create dimensional plots. These parameters are in the appendix of the paper.

    Figure production for paper. Imports text files and uses parameters from param.m. Each figure is produced by excecuting a section of code.

    Contains the data exported from COMSOL required to run 4).
