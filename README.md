# steadystate
Matlab-based semi-automated routine to detect the optimum combination of first new steady-state friction point following a friction velocity step and the associated slip window to accurately operate a linear detrend starting the new steady-state point in order to determine accurate rate- and state- friction parameter values from inverse modelling.

Branch main: contains 1) the code required to slice the friction tests into its single friction velocity steps (slicing_velsteps.m); 2) the main code to determine steady-state friction after the velocity steps (steadystate.m); 3) the subroutine of the the main code steadystate.m (getpoints_velstep.m)

Branch Demo-tests: contains a .txt file to test the code slicing_velstep.m (TR2_exp013_sepiolite_to-be-sliced.txt) and n.4 single velocity steps that can be run directly with the code steady-state.m (i.e., velstep1,2,3_TR2_exp013_sepiolite.txt and velstep7_BRAVA_expb819_basalt_Giacomel_et_al_2021.txt 
