function [Time,Disp,Mu,sneff,LS2,steady_state] = steadystate(min_LS2,max_LS2,delta_LS2,deltaslip); %semicomma to suppress outputs

%% FUNCTION TO DETERMINE THE FIRST STEADY-STATE POINT FOLLOWING A FRICTION VELOCITY STEP and the associated regression window length starting that point
%% Program description
%  The approach followed in this routine for determining new steady-state conditions, is the systematic estimation of slopes after the velocity
%  jump (S2) via linear regression, until S2 approaches the slope before the velocity step (S1) such that S2 ~ S1.
%  Slope S2 is systematically estimated within an evaluation slip window (LS2), which is displaced by the amount deltaslip towards larger displacements, 
%  until the slope S2 reaches or falls below the threshold of 5e-7 microns^-1. 
%  The displacement corresponding to the first point satifying this condition, associated with a given LS2 is called "first new steady-state point". 

%  This routine allows to evaluate the first new steady-state point using multiple slip window lengths from 'min_LS2' to
%  'max_LS2' with step 'delta_LS2', and suggests the user an optimum combination of first new steady-state point and slip window length LS2. 
%  This pair of values should be used for the preliminary linear detrend operation preceding the inverse modelling required for getting the rate- and state- friction (RSF) parameters.

%  In velocity steps where the slope S1 is larger than 5e-7 microns^-1, 
%  S1 is removed internally in the program to be able to work with the same convergence criterion (i.e., |S2|<=5e-7 microns^-1) in all the cases. 

%  The routine works well provided that steady-state conditions have been reached before the velocity step.

%% Input parameters
%  min_LS2: minimum size of the slip window (microns) within which slope S2 has been worked out after the velocity step; 
%           if omitted, by default min_LS2 = 50 microns.
%  max_LS2: maximum size of the slip window (microns) within which slope S2 has been worked out after the velocity step;
%           max_LS2 by default is total window length (TWL)-dependent, where TWL is the slip interval from the velocity step+200 microns to the last selected point of the datafile 
%           if max_LS2 is omitted, by default: if TWL <= 1000 microns, max_LS2 = 50% of the approximated (to the multiples of a hundred) total window length (TWL). 
%                                              if TWL > 1000 microns,  max_LS2 = 500 microns                                        
%  delta_LS2: increment in size of the slip window (microns) from one steady state-analysis to the next one; 
%             if omitted, by default delta_LS2 = 50 microns. 
%  deltaslip: amount indicating how much the slip window is displaced from one slope S2 estimation to the next one within a steady-state analysis; 
%             if entered by the user, is in microns. For example, a deltaslip = 10 means that the user is calculating the slope S2 within a given slip window length LS2 and moving such a window every 10 microns;
%             if omitted, by default deltaslip = equivalent slip from one datapoint to the next one after 200 microns past the velocity step

%% Output parameters
%  Time, Disp, Mu, sneff: time,displacement,friction, and effective normal stress data, respectively, to be entered in the rate- and state- friction analysis that follows this one 
%  (data have already been already windowed, from points 1 to 3 selected by the user using the function getpoints_velstep.m, with index specified by the DataIndex array) 
%  
%  LS2: array of slip windows lengths (microns) used to estimate S2, that spans from min_LS2 to max_LS2 with increment of delta_LS2 from one analysis to the next one
%  steady_state: array of 1st new steady-state points (microns) associated with a given slip window length LS2 used for the steady-state analysis 

%% How to run the code:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE USE OF DEFAULT SETTINGS IS HIGHLY RECOMMENDED %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - To use the default input parameters without displaying the output parameters in the workspace: 
%   press the MATLAB Run button or alternatively type: 'steadystate;' in the MATLAB Command Window
% - To use the default input parameters and record the outputs in the workspace (suggested option when using the MATLAB program RSFit3000 to determine the RSF parameters https://github.com/rmskarbek/RSFit3000):
%   [Time,Disp,Mu,sneff,LS2,steady_state]= steadystate;

% - To enter customized parameters: replace values for the following input paramenters: steadystate(min_LS2,max_LS2,delta_LS2,deltaslip)
%   After entering the values, delta_LS2 will be rounded to the multiple of tens and in turn, min_LS2 and max_LS2 will be approximated accordingly. 

% example: change all default parameters:
% [Time,Disp,Mu,sneff,LS2,steady_state] = steadystate(20,360,20,10); 
% where min_LS2=20, max_LS2=360, delta_LS2=20 and deltaslip=10.

% When both customised and default input arguments coexist, use [] or do not input the arguments for the default settings:

% example1: keep deltaslip as the default and change the others:
% [Time,Disp,Mu,sneff,LS2,steady_state] = findsteadystate(10,100,20,[]) 
% or equivalently,
% [Time,Disp,Mu,sneff,LS2,steady_state] = findsteadystate(10,100,20)
% example2: change only deltaslip and keep the others as default: [Time,Disp,Mu,sneff,LS2,steady_state] = steadystate([],[],[],10);
% example3: change min_LS2 and max_LS2 and keep the others as default: [Time,Disp,Mu,sneff,LS2,steady_state] = steadystate(10,100);

%% Upload experimental parameters within the vel step, i.e.,DataIndex, Time, Displacement, Friction and normal stress 
a = input('Have you already sliced your experimental data into single velocity steps(y/n)? \n','s');
if a=='n'
   error('Slice your experiments e.g. with slicing_velsteps.m before running the program')
end
if a=='y' % import single velocity step
   [DataIndex,Time,Disp,Mu,sneff]=getpoints_velstep; %call to getpoints_velsteps.m
end

%% Default checking
%Input parameters
close all; clc;
disp_pt3 = Disp(DataIndex(2))+200; % start displacement for the slope S2 estimation: point 2 + 200 microns
ind_pt3=find(Disp>=disp_pt3); ind_pt3= ind_pt3(1); % find DataIndex associated with slip = slip(DataIndex(2))+ 200 microns

TWL = Disp(DataIndex(3))-Disp(ind_pt3); %total window length available for the slope S2 estimations
decimal = -2; roundTWL = round(10^decimal*TWL)/(10^decimal); %approximate TWL to the nearest multiple of a hundred 

if nargin < 4 || isempty(deltaslip)
    n=1; % n = 1 means that slope S2 is calculated every single data point  
%     checkslipint = Disp(ind_pt3+n)-Disp(ind_pt3); %check delta displacement right after the velocity step
else
    avgdeltadisp = abs(mean(diff(Disp(ind_pt3:DataIndex(3))))); %average slip distance from one data point to the next one starting from the velocity step + 200 microns until the end of the datafile
    if (deltaslip<avgdeltadisp)  %if the deltaslip inputted by the user is lower than 
                                 %the actual average delta displacement 'avgdeltaslip', the program sets by default 
                                 %the deltaslip interval equal to the minimum slip distance between two points and warns the user about the changes made
        n=1;
        warning ('input deltaslip is smaller than the slip distance between two consecutive points=');
        text = ['deltaslip has been set by default equal to the minimum delta slip = ',num2str(avgdeltadisp),' ',char(956),'m'];
        disp(text)
%         checkslipint = Disp(ind_pt3+n)-Disp(ind_pt3); %check the accuracy of the conversion in data points of deltaslip
    else
        n=round(deltaslip/avgdeltadisp); %conversion of deltaslip (microns) in equivalent datapoints
%         checkslipint = Disp(ind_pt3+n)-Disp(ind_pt3); %check the accuracy of the conversion in data points of deltaslip
    end
end 
if nargin < 3 || isempty(delta_LS2)
   delta_LS2 = 50; % microns    
else
     decimal1= -1; delta_LS2 = round(10^decimal1*delta_LS2)/10^decimal1; %round up delta_LS2 to the closest multiple of ten when entered by the user
end
if ((nargin < 2 || isempty(max_LS2))&&(TWL<=1000))
    max_LS2 = 0.5*roundTWL; %when TWL <= 1000 microns, by default max_LS2 = 0.5*TWL rounded up to the multiple of hundreds microns  
elseif ((nargin < 2 || isempty(max_LS2))&&(TWL>1000))
    max_LS2 = 500; %when TWL > 1000 microns, by default  max_LS2 = 500 microns
else
    ratiomaxLS2 = max_LS2/delta_LS2; decimal2=0; roundmaxratio = round(10^decimal2*ratiomaxLS2)/10^decimal2; 
    max_LS2 = delta_LS2*roundmaxratio; %round intervals entered by the user in order to have multiples of delta_LS2 for max_LS2
end
if nargin < 1 || isempty(min_LS2)
   min_LS2 = 50; % microns
else
    ratiominLS2 = min_LS2/delta_LS2; decimal2=0;roundminratio = round(10^decimal2*ratiominLS2)/10^decimal2; 
    min_LS2 = delta_LS2*roundminratio; %round intervals entered by the user so that to have multiples of delta_LS2 for min_LS2
end

txtmin_LS2=['min_LS2 = ',num2str(min_LS2),' ',char(956),'m'];
disp(txtmin_LS2);
txtmax_LS2=['max_LS2 = ',num2str(max_LS2),' ',char(956),'m'];
disp(txtmax_LS2);
txtdelta_LS2=['delta_LS2 = ',num2str(delta_LS2),' ',char(956),'m'];
disp(txtdelta_LS2);
fprintf('\n');

% when max_LS2 is entered by the user, it displays in the command window the most recommended value based on the inputted max_LS2 and TWL  
if (max_LS2>0.5*roundTWL&&TWL<=1000) 
    suggested_maxLS2 = 0.5*roundTWL;
    suggestion = ['suggested max_LS2 = ',num2str(suggested_maxLS2),' ',char(956),'m'];
    disp(suggestion);
    warning ('max_LS2 is larger than 50% of the approximated total slip window length: it is advised to rerun the routine by inputting max_LS2 as suggested above');
elseif (max_LS2>0.5*roundTWL&&TWL>1000)
    suggested_maxLS2 = 500;
    suggestion = ['suggested max_LS2 = ',num2str(suggested_maxLS2),' ',char(956),'m'];
    disp(suggestion);
    warning ('max_LS2 is larger than 50% of the approximated total slip window length, which is larger than 1000 microns: it is advised to rerun the routine by inputting max_LS2 as suggested above');
end

detr = (min_LS2:delta_LS2:max_LS2)'; %intervals of slip windows LS2 used for estimating S2 slopes throughout the multiple steady-state analyses
if length(detr)<3
    error ('Please set the number of slip window length (LS2) intervals to a minumum of 3 values (e.g., min_LS2 = 100, max_LS2 = 200, ,delta_LS2 = 50)');
end

%% Plot velocity step from point 1 to 3
figure(1)
plot(Disp(DataIndex(1):DataIndex(3)),Mu(DataIndex(1):DataIndex(3)),'-k','Linewidth',1); %plot the friction vs. displacement event from the selected point 1 to point 3
ax=gca; ax.XAxis.Exponent = 0;
xlabel ('Displacement (\mum)');
ylabel ('Friction coefficient');hold on;

%% PHASE 1: Slope S1 estimation prior to vel step between DataIndex 1 and 2 (i.e., from the starting point to the one where the velocity jump occurs - all the points here are at steady state) 
format long

tdata_bvs = Time(DataIndex(1):DataIndex(2)); %time (s) in seconds before the velocity step
xdata_bvs = Disp(DataIndex(1):DataIndex(2)); %displacement (microns) before the velocity step
ydata_bvs = Mu(DataIndex(1):DataIndex(2));  %friction coefficient at steady-state conditions before the velocity step
 
linreg1 = fitlm(xdata_bvs,ydata_bvs); %linear regression before the velocity step
noise_bvs = linreg1.RMSE; %square root of the mean squared error of the linear regression (root mean square error) - proxy for the noise level of the signal before the velocity step
slope1= linreg1.Coefficients.Estimate(2); %slope before the velecity step S1

dispint = Disp(DataIndex(2))-Disp(DataIndex(1)); %displacement interval before the velocity step where slope S1 is estimated

deltat_bvs = mean(diff(tdata_bvs)); %average delta time (s) before the velocity step 
freq_bvs = 1/deltat_bvs; %average sampling frequency (Hz) before the velocity step 
deltadisp_bvs = abs(mean(diff(Disp(DataIndex(1):DataIndex(2))))); %average delta displacement (microns) before the velocity step
vel_bvs = deltadisp_bvs/deltat_bvs; %average slip velocity (microns/sec) before the velocity step 
Npointsperdisp_bvs = freq_bvs/vel_bvs; %number of points per unit displacement (microns^-1) before the velocity step = freq/vel = 1/deltadisp 

% Check slope S1 accuracy
% Errors based on slope S1 estimations from a friction database containing  synthetic velocity steps from V1 = 0.3 to V2 = 3 microns/s 
% with no overall slip  dependencies in friction (S1~0), generated with the following input parameters: a = 0.012; b1 = 0.010; b2 = 0.007; DRS_1 = 20 microns; 
% DRS_2 = 500 microns; k'= 0.01 microns^-1. Each velocity step in the database, differs in the noise level  
% SD (~RMSE in zero-slope friction data) = 0.00025-0.000375-0.0005-0.00075-0.001, and the sampling frequency f = 0.1-1-10-100-1000 Hz. 

% When in the database combinations of noise and f result in estimations of slope S1 |S1|>5x10^-7 microns^-1, S1 determinations are regarded as erratic.  
% Due to this, the routine sends an error when comparable levels of RMSE and number of no. points/unit slip (=f/V) are reached in the experimental data as likely detrimental to slope S1 estimations.
% In this routine, a 10 % tolerance has been applied to the noise level, the no. points/unit slip, and the slip window length (LS1) used for S1 estimations
% contained in the database.

if  dispint<45
    disp('flag 5')
    warn=['slope1 = ',num2str(slope1),' 1/',char(956),'m'];
    disp(warn)
    fprintf('\n');
    error('The slope measurement before the velocity step might not be accurate, hence the whole steady-state assessment: it is advised to choose either a larger slip interval for the linear regression before the vel step (LS1) or run a test with higher sampling frequency/larger step length');
end
if ((dispint<67.5&&dispint>=45)&&(noise_bvs<=0.000825)&&(Npointsperdisp_bvs<30))||((dispint<67.5&&dispint<=45)&&(noise_bvs>0.000825)&&(Npointsperdisp_bvs<3000))
    disp('flag 4')
    warn=['slope1 = ',num2str(slope1),' 1/',char(956),'m'];
    disp(warn)
    fprintf('\n');
    error('The slope measurement before the velocity step might not be accurate, hence the whole steady-state assessment: it is advised to choose either a larger slip interval for the linear regression before the vel step (LS1) or run a test with higher sampling frequency/larger step length');
end 
if ((dispint<90&&dispint>=67.5)&&(noise_bvs<=0.0004125)&&(Npointsperdisp_bvs<3))||((dispint<90&&dispint>=67.5)&&(noise_bvs<=0.000825&&noise_bvs>0.0004125)&&(Npointsperdisp_bvs<30))||((dispint<90&&dispint>=67.5)&&(noise_bvs>0.000825)&&(Npointsperdisp_bvs<300))
    disp('flag 3')
    warn=['slope1 = ',num2str(slope1),' 1/',char(956),'m'];
    disp(warn)
    fprintf('\n');
    error('The slope measurement before the velocity step might not be accurate, hence the whole steady-state assessment: it is advised to choose either a larger slip interval for the linear regression before the vel step (LS1) or run a test with higher sampling frequency/larger step length');
end 
if ((dispint<112.5&&dispint>=90)&&(noise_bvs<=0.0004125&&noise_bvs>0.000275)&&(Npointsperdisp_bvs<3))||((dispint<112.5&&dispint>=90)&&(noise_bvs>0.0004125)&&(Npointsperdisp_bvs<30))
    disp('flag 2')
    warn=['slope1 = ',num2str(slope1),' 1/',char(956),'m'];
    disp(warn)
    fprintf('\n');
    error('The slope measurement before the velocity step might not be accurate, hence the whole steady-state assessment: it is advised to choose either a larger slip interval for the linear regression before the vel step (LS1) or run a test with higher sampling frequency/larger step length');
end 
if (dispint>=112.5&&noise_bvs>0.0004125&&Npointsperdisp_bvs<3)
    disp('flag 1')
    warn=['slope1 = ',num2str(slope1),' 1/',char(956),'m'];
    disp(warn)
    fprintf('\n');
    error('The slope measurement before the velocity step might not be accurate, hence the whole steady-state assessment: it is advised to choose either a larger slip interval for the linear regression before the vel step (LS1) or run a test with higher sampling frequency/larger step length');
end 

%Refinement of slope S1: removal of outliers using the studentized residuals
stud_res = linreg1.Residuals.Studentized; %determination of the studentized residuals to identify outliers and exclude them from the slope analysis before the vel step
I_studres3 = abs(stud_res) > 3; %outliers: all points outside +-3 sigma
outliers_studres3 = excludedata(xdata_bvs,ydata_bvs,'indices',I_studres3); %99.7% of the data are within 3 standard deviations of the mean
linreg1_filt = fitlm(xdata_bvs,ydata_bvs,'Exclude',outliers_studres3); %linear regression before the velocity step following the removal of outliers 
slope1_filt= linreg1_filt.Coefficients.Estimate(2); %estimation of slope S1 before the velecity step after excluding the outliers

%Preprocessing of the velocity step for consistent steady-state analysis independently of overall slip dependencies in friction (done only internally in the script) 
if  abs(slope1_filt)>5e-7 % if |S1|>5e-7 microns^-1
    Mudetr = Mu+slope1_filt.*(Disp(DataIndex(2))-Disp); %pre-processing the data: removing the linear trend before the velocity step using the slope S1 
else
    Mudetr = Mu; % if |S1|<=5e-7 microns^-1 no need to linear detrend the data
end

%% PHASE 2: Calculate the slopes after the velocity step and find the steady-state points at a given slip window length LS2

tdata_avs = Time(ind_pt3:DataIndex(3)); %time (s) after the velocity step
deltat_avs = mean(diff(tdata_avs)); %average delta time (s) after the velocity step  
freq_avs = 1/deltat_avs; %average sampling frequency (Hz) after the velocity step  
deltadisp_avs = abs(mean(diff(Disp(ind_pt3:DataIndex(3))))); %delta displacement (microns) after the velocity step
vel_avs = deltadisp_avs/deltat_avs; %slip velocity (microns/s) after the velocity step 
Npointsperdisp_avs = freq_avs/vel_avs; %number of points per unit displacement (microns^-1) = 1/deltadisp; after the velocity step

warnings = 0; %initialize counter for warnings that arise in the case that the first steady-state point found by the script is the last available data point in the dataset 

for j = 1:length(detr) %reiterate the steady-state analysis for every slip window length from min_LS2 to max_LS2 with increment delta_LS2 from one analysis to the next one
    txt = ['Linear regression window size = ',num2str(detr(j)),' ',char(956),'m'];
    disp(txt)
    lastptevalwin = Disp(DataIndex(3))-detr(j); Nlastptevalwin = find(Disp>lastptevalwin,1)-1; Ndetr_end= length(Nlastptevalwin:DataIndex(3)); %last evaluation slip window for slope S2 in equivalent data points
    lastevalptdetr = Disp(ind_pt3)+detr(j); Nlastevalptdetr = find(Disp>lastevalptdetr,1)-1; Ndetr = length(ind_pt3:Nlastevalptdetr); %evaluation slip windows throughout the analysis in equivalent data points 
    %we made the differentiation between Ndetr_end and Ndetr since it is possible that the delta displacement might vary along the velocity step due to slight changes in velocity, 
    %hence, the number of data points to get the same displacement may change throughout the test.
    
    evalwin = (ind_pt3:(DataIndex(3)- Ndetr_end))'; %slip window length in equivalent data points for the calculation of the slopes S2 after the velocity step  
    lastind = 1+((fix((length(evalwin)-1)/n))*n); %generalization of the last datapoint (in data indexes) within evalwin for the slope S2 analysis

    %To quantify the mean noise (RMSE) after the velocity step (through a ~200 microns slip window towards the end of the velocity step)
    startlinregr4noise = DataIndex(3)-round(Npointsperdisp_avs*250); endlinregr4noise = round(DataIndex(3)-Npointsperdisp_avs*50); %intervals of the linear regression model used to quantify the noise after the velocity step
    linregress4noise = fitlm(Disp(startlinregr4noise:endlinregr4noise),Mudetr(startlinregr4noise:endlinregr4noise));
    noise_avs = linregress4noise.RMSE; %square root of the mean squared error of the linear regression (root mean square error) - proxy for the noise level of the signal after the velocity step

    lindetr2 = polyfit(Disp(ind_pt3:ind_pt3+Ndetr),Mudetr(ind_pt3:ind_pt3+Ndetr),1); %1st linear regression after velocity step + 200 microns
    slope2 = lindetr2(1); %slope S2 associated with the first regression line after velocity step + 200 microns
    
    steadystate(j)=Disp(evalwin(1)); % initialize steady-state: the first steady-state candidate (microns) is the first point of the first regression line 
                                     % used to estimate slope S2 after the velocity step (i.e., DataIndex2 + 200 microns) 
    
    i=1; % initialize counter for slope S2 estimation
       
    while (abs(slope2)>5e-7&&i<lastind) %iterative loop to find steady state continues until |S2|<=5e-7 microns^-1 or until the last available datapoint is reached
          i=i+n; %update index i by increasing it of the amount n during each cycle of the iterative loop          
          lastevalptdetr = Disp(evalwin(i))+detr(j); Nlastevalptdetr = find(Disp>lastevalptdetr,1)-1; Ndetr= length(evalwin(i):Nlastevalptdetr);  
          lindetr2 = polyfit (Disp(evalwin(i):(evalwin(i)+Ndetr)),Mudetr(evalwin(i):(evalwin(i)+Ndetr)),1); slope2 = lindetr2(1); 
          % update slope S2 calculation as you move towards larger displacements along the velocity step
          steadystate(j) = Disp(evalwin(i)); %update steady-state candidate point (microns)   
          %slopetwodetr(j)=slope2(end); %to store the last slope S2 relative to the steady-state analysis               
    end   
    if   (steadystate(j)<Disp(evalwin(lastind)))
          results = ['Steady state = ',num2str(steadystate(j)),' ',char(956),'m'];
          disp(results); %display results in the MATLAB Command Window if steady-state condition is met before the last available point of analysis
    else
        warnings = warnings+1; %if steady-state condition is met upon the last available point of analysis, update the warming counter and display a warning on this
        warning ('Steady-state condition never met or corresponds to the last point available for the analysis and thus might not be accurate');
    end   
    
    fig=figure(101); %retrieve from the figure the x and y data (i.e., LS2 and first new steady-state)
    h(j)=plot(detr(j),steadystate(j));
    set(fig,'position',[-500,-500,1,1]); 
    set(fig,'visible','off'); % in order not to visualize the plot
    x(j) = get(h(j),'XData');LS2=x';
    y(j) = get(h(j),'YData');steady_state=y'; 

    figure(99) %plot outputs of the analyses: first new steady state vs. slip window length used to estimate slope S2 (LS2)
    plot(LS2,steady_state,'-ko','markersize',8,'markerfacecolor','k');
    ax=gca; ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0;
    xlabel ('Linear regression window length, L_{S2} (\mum)');
    ylabel ('First new steady-state point (\mum)');
    xticks([detr])
    xlim([detr(1)-delta_LS2 detr(end)+delta_LS2]);
    ylim([Disp(ind_pt3) Disp(DataIndex(3))])
    
    if (warnings >1) % if the warning relative to the last datapoint available appears more than once, close the program and display the error with the suggestions
          clc;close all
          txtmin_LS2=['min_LS2 = ',num2str(min_LS2),' ',char(956),'m'];
          disp(txtmin_LS2);
          txtmax_LS2=['max_LS2 = ',num2str(max_LS2),' ',char(956),'m'];
          disp(txtmax_LS2);
          txtdelta_LS2=['delta_LS2 = ',num2str(delta_LS2),' ',char(956),'m'];
          disp(txtdelta_LS2);
          fprintf('\n');
          if ((max_LS2>0.5*roundTWL))
             disp(suggestion);
          end
          error ('Steady-state condition never met or corresponds to the last point available for the analysis and thus might not be accurate. Run a new friction velocity test with longer step length.')              
    end
end

%display outputs in the command window
fprintf('\n');
txt1 = ['Average noise level after the velocity step = ',num2str(noise_avs)];
txt2 = ['Number of points per unit displacement after the velocity step = ',num2str(Npointsperdisp_avs)]; %displays in the command window the noise level and freq/vel
disp(txt1)
disp(txt2)

%% Among the first new steady-state points determined for each given LS2, suggest the optimum pair of 1st new steady state-LS2 for linear detrending in RSF software (steadybest in the script)
steadysort = sort(steady_state); %Sort results in ascending order
steadybest = []; %initialize steadybest

diffsteady = abs(diff(steady_state)); %differential of steady-state points
diffLS2 = abs(diff(LS2)); %differential of linear regression window lengths used for S2 slope estimations

% Find the first new steady-state points that don't change significantly with the linear regression window length LS2 
% and among the largest values, find the one associated with the smallest LS2. 
% If this point exists, this is the called the optimum first new steady-state point coupled with LS2.

c = find(abs(diffsteady./diffLS2)<=1.5); % find optimum 1st steady state point: 
% STEP 1)calculate the ratio between the differential of first steady-state points and the differential relative to LS2; 
% 1.5 is a threshold value based on experience on several velocity steps when steady state was evaluated using deltaLS2 = 50 microns 
d = find(steady_state(c)>=steadysort(length(steadysort)-2)); %STEP 2) find the first 3 largest steady-state points satisfying the above condition

if (length(d)>0)
    steadybest = steady_state(c(d(1))); % if there is a candidate, the optimum first steady-state point is the first one on the list, i.e., the value associated with the smallest LS2 among the steady-state points
    fprintf('\n'); %display combination of optimum first new steady-state point and slip window length LS2 suggested for performing the linear detrend during the subsequent RSF analysis
    txt3 = ['Optimum combination of 1st steady-state point = ',num2str(steadybest),' ',char(956),'m'];
    txt4 = ['and the associated window size to linear detrend = ',num2str(LS2(c(d(1)))),' ',char(956),'m'];
    disp(txt3);
    disp(txt4);
    fprintf('\n');
    
    figure(); %plot outputs from the routine with highlighted the optimum pair of first new steady state - LS2
    plot(LS2(c(d(1))),steady_state(c(d(1))),'-o','markersize',15,'MarkerFaceColor','r','markeredgecolor','r');hold on;
    plot(LS2(c(d(1))),steady_state(c(d(1))),'o','markersize',12,'MarkerFaceColor','k');hold on;
    plot(LS2,steady_state,'-o','markersize',8,'MarkerFaceColor','g','linewidth',1)
    axes = gca; axes.XAxis.Exponent = 0;  axes.YAxis.Exponent = 0; axes.ColorOrder = [0 0 0];
    xlabel ('Linear regression window length, L_{S2} (\mum)');
    ylabel ('First new steady-state point (\mum)');
    xticks([detr])
    xlim([detr(1)-delta_LS2 detr(end)+delta_LS2]);
    ylim([Disp(ind_pt3) Disp(DataIndex(3))])
    close(99);
    
    figure(1); %alongside the 1st new steady state vs LS2 graph, display the original velocity step with superimposed the optimum 1st new steady state and the associated detrend line LS2 
    x1=xline (steady_state(c(d(1))),'--r',{'STEADY-STATE'});x1.LabelHorizontalAlignment ='left';hold on;
    x2=xline (steady_state(c(d(1))),'--r',{'Linear', 'detrend', 'here'});hold on;
    x2.LabelOrientation ='horizontal';x2.LabelVerticalAlignment ='top';x2.LabelHorizontalAlignment ='right';hold on;
    xline(steady_state(c(d(1)))+LS2(c(d(1))),'--r');
    
else
   close all;clc;
   txtmin_LS2=['min_LS2 = ',num2str(min_LS2),' ',char(956),'m'];
   disp(txtmin_LS2);
   txtmax_LS2=['max_LS2 = ',num2str(max_LS2),' ',char(956),'m'];
   disp(txtmax_LS2);
   txtdelta_LS2=['delta_LS2 = ',num2str(delta_LS2),' ',char(956),'m'];
   disp(txtdelta_LS2);
   fprintf('\n');
   error('No optimum new steady-state points. Repeat the analysis by increasing LS2_max of a displacement amount equal to delta_LS2: type "steadystatetry(min_LS2,max_LS2)" and enter new numeric value for max_LS2. Alternatively, run a new velocity step with higher frequency/slip velocity ratio');
end    

%% Save outputs in .xlsx format
a = input('Save outputs(y/n)? \n','s');
    
%save in the workspace
Time = Time(DataIndex(1):DataIndex(3)); %Slicing of time, displacement, friction, and effective normal stress data using the first and the last selected points in order 
Disp = Disp(DataIndex(1):DataIndex(3)); %to input them in the MATLAB software used to retrieve the modelled rate- and state-friction parameter values (https://github.com/rmskarbek/RSFit3000)
Mu = Mu(DataIndex(1):DataIndex(3));
sneff = sneff(DataIndex(1):DataIndex(3));
    
if a=='n'
end
if a=='y'
    %save main outputs (first new steady-state points vs. linear regression slip window LS2) in the current folder in Excel spreadsheets
    filename=input('Save outputs: please input file name.xlsx\n','s');
    detrendtxt = ['LS2(',char(956),'m)']; steadytxt = ['1st new steady-state(',char(956),'m)'];
    hpluto = {detrendtxt   steadytxt};
    cpluto = [LS2 steady_state];
    cpluto_opt = [LS2(c(d(1))),steady_state(c(d(1)))];
    Mpluto = [hpluto;num2cell(cpluto)];
    Mpluto_opt = [hpluto;num2cell(cpluto_opt)];
    if (length(d)>0)
        xlswrite(filename,Mpluto,'all outputs') % save in the first spreadsheet of the Excel file the compilation of outputs
        xlswrite(filename,Mpluto_opt,'optimum output') % save in the second spreadsheet of the Excel file the optimum output
    else
        xlswrite(filename,Mpluto,'all outputs') 
    end 
end

%% Overwrite original .txt file with data from selected point 1 to 3

s = input('Save datafile in txt format from selected points 1 to 3 (y/n)? \n','s');

if s=='n'
end
if s=='y'
    name=input('Enter name of the velocity step and experiment.txt\n','s');
    fileID = fopen(name,'w');
    fprintf(fileID,'%s\t%s\t%s\t%s\n','Time(s)','Displ(microns)','Friction','Effective_Normal_stress(MPa)');             
    %headers of each .txt file
    fprintf(fileID,'%f\t%f\t%f\t%f\n',[Time Disp Mu sneff]');
    fclose(fileID); 
end
end