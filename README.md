Material for paper: 
"A continuum model for the elongation and orientation of Von
Willebrand Factor with applications in arterial flow"

Link to archive article: 

Author of files in this repository: E. Yeo, email edwina.yeo.14@ucl.ac.uk

Contents:

1. Python files (requires fenics2019.2.0.dev0 which is imported through dolfin):

    main.py
    main file runs all simulations for data presented in paper. Each figure is produced by a different 'study' appropriately labeled. Data for each figure is plotted in 	the results section. Data was then exported to MATLAB to provide uniform formatting and visualisation.

    Numerical verification file to suport main file. Model unchanged, mesh convergence is tested through simulations with increasingly fine meshes. Results sensitivity to 	the approximate parameter 'delta' is tested.

    Contains parameters used to create dimensional plots. These parameters are in the appendix of the paper.


2. GMSH files (requires gmsh-3.0.6):
	GMSH files are generated and meshed in python script main.py


3. MATLAB Files (requires MATLAB lisence) 3) param.m 4) plot_main_text.m 5) Data
	Figure production for paper. Imports text files and uses parameters from param.m. Each figure is produced by excecuting a section of code. Figure 4 is augmented 	via fene.svg using Inkscape.

4. Raw Data exported python stored in Data folder

5. Inkscape SVG files for diagrams and Figure 4 augmentation.




    
