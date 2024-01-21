# steadystate
In this GitHub repository, you will find the MATLAB-based automated routine that enables the determination of the optimum combination of the first new steady-state friction point following a velocity step and the associated slip window length required to accurately operate a linear detrend starting from that point. This preliminary operation is key especially for materials characterized by long characteristic slip distances D_{RS} (i.e., large frictional transients following the velocity steps), to determine accurate rate- and state- friction parameter values via inverse modelling.

CONTENT OF THE REPOSITORY:
1) Branch "main":

It consists of two main files:
the main script steadystate.m and its subscript getpoints_velstep.m.
In order to run the routine, both files must be on the same MATLAB path. 
The routine can be run using the default settings (recommended option) by typing steadystate in the MATLAB command window or by pressing Run in the EDITOR window of the steadystate.m file opened in MATLAB.
The input and output parameters of steadystate.m are fully described in lines 20 to 41 of the m-file.
For customised input parameters and relative examples we refer to lines 96 to 130 of steadystate.m. 

Since steadystate can process only one velocity step at a time, it requires pre slicing of the friction datafile, which can be done using the script slicing_velsteps.m, also attached in the main branch. 
This is a straightforward, standalone script that can be run like any other m-file.
Each sliced velocity step is graphed at the end of the slicing (friction vs. time and friction vs. displacement curves) and stored as Time, Displacement, Friction, and Effective normal stress arrays in a .txt file automatically generated in the selected MATLAB path.

2) Branch Demo-velocity-stepping-tests: contains a .txt file to test the code slicing_velstep.m (i.e., TR2_exp013_sepiolite_to-be-sliced.txt) and n.3 single velocity steps that can be run directly with the code steadystate.m, of which n.1 is synthetic (i.e., synthetic_D_{RS2}_500microns.txt) and n.2 are from experimental datasets of different rock deformation laboratories (i.e., velstep3_TR2_exp013_sepiolite.txt; velstep7_BRAVA_expb819_basalt.txt). 


# Citing this work
If you find these codes helpful and use them for your projects, please consider citing the following paper describing the underlying steady state method:

Giacomel et al., 2024 - GSA, Geosphere 
"steadystate: A MATLAB-based routine for determining steady-state friction conditions in the framework of rate- and state- friction analysis."

