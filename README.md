# steadystate
In this GitHub repository, you will find the MATLAB-based automated routine that enables the determination of the optimum combination of the first new steady-state friction point following a velocity step and the associated slip window length required to accurately operate a linear detrend starting from that point. This preliminary operation is required to determine accurate rate- and state- friction parameter values from inverse modelling.

Branch main: 
steadystate consists of two main files:
the main script steadystate.m and its subscript getpoints_velstep.m.
In order to run the routine, both files must be on the same MATLAB path. 
The routine can be run using the default settings by typing steadystate in the MATLAB command window or by pressing Run in the EDITOR subfolder of the steadystate.m file opened in MATLAB.

The program along with its input and output parameters are fully described in lines 4 to 38 of the steadystate.m file.
For customised input parameters and relative examples we refer to lines 44 to 63 of steadystate.m. 

Since steadystate can process one velocity step at a time, it requires pre slicing of the friction datafile, which can be done using the script slicing_velsteps.m, also attached in the main branch. 
This is a straightforward script that can be run like any other .m file.
Each sliced velocity step is graphed at the end of the slicing (friction vs. time and friction vs. displacement curves) and stored as Time, Displacement, Friction, and Effective normal stress arrays in a .txt file in the selected MATLAB path.

Branch Demo-tests: contains a .txt file to test the code slicing_velstep.m (TR2_exp013_sepiolite_to-be-sliced.txt) and n.5 single velocity steps that can be run directly with the code steadystate.m (i.e., velstep1,2,3_TR2_exp013_sepiolite.txt and velstep4,7_BRAVA_expb819_basalt_Giacomel_et_al_2021.txt).

# Citing this work
If you find these codes helpful and use them for your projects, please consider citing the following paper describing the underlying steady state method:

Giacomel et al. steadystate: A MATLAB-based routine for determining steady-state friction conditions in the framework of rate- and state- friction analysis.

Currently under review in GSA, Geosphere. 
