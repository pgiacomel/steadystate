function [Time,Disp,Mu,sn,LS2,steady_state] = steadystate(deltaslip,min_LS2,delta_LS2,max_LS2); %semicomma to suppress outputs

%% FUNCTION TO DETERMINE THE FIRST STEADY-STATE POINT FOLLOWING A FRICTION VELOCITY STEP and the associated regression window length starting that point
%% Program description
%  The approach followed in this program for determining steady-state condition at the new slip velocity of a friction velocity step test (= new steady-state), 
%  is that the slope estimated from linear regression after the velocity step (S2) approaches the slope before the velocity step (S1), i.e., S1 ~ S2.
%  Slope (S2) is systematically estimated within an evaluation slip window (deltaslip), which is displaced along the velocity steps towards larger displacements, 
%  until the slope S2 reaches or falls below the threshold of 5e-7 microns^-1. 
%  The displacement corresponding to the first point satifying this condition is called "first new steady-state point". 
%  When the slope before the velocity step S1 is larger than 5e-7 microns^-1, 
%  S1 is removed internally in the program to be able to work with the same convergence criterion (i.e., |S2|<5e-7 microns^-1) for steady-state 
%  in all the cases, also those in which slip dependency in friction before the velocity step is significantly different than 0.

%  This program allows to evaluate the first new steady-state point using multiple moving slip window lengths (from 'min_LS2' to
%  'max_LS2' with step 'delta_LS2'), and suggests the user the optimum combination of first new steady-state point and slip window length that should be used 
%  for the preliminary detrend operations preceding the inverse modelling required for getting the rate- and state- friction (RSF) parameters.
  
%  The program works well assuming that steady-state conditions have been reached before the velocity step, which is one of the essential conditions for
%  determining good modeled RSF parameter values.

%% Input parameters
%  deltaslip: amount indicating how much the slip window is displaced from one slope calculation after the vel step to the other; 
%             for example, a deltaslip = 50 means that the user is calculating the slope within a given moving slip window length LS2 and moving such window every 50 microns;
%             if omitted, by default deltaslip = 1 data points (which can differ from 1 micron) 
%  min_LS2: minimum size of the moving slip window (microns) within which slope has been worked out after the velocity step 
%           if omitted, by default min_LS2 = 50 microns
%  delta_LS2: increment in the size of the moving slip window (microns); 
%           if omitted, by default delta_LS2 = 50 microns 
%  max_LS2: maximum size of the moving slip window (microns) within which slope has been worked out after the velocity step
%           if omitted, by default max_LS2 = 500 microns if the velocity step length is >= 1000 microns
%                                                      = 300 microns if the velocity step length is < 1000 microns 

%% Output parameters
%  Time Disp Mu sn: respectively: time,displacement,friction and normal stress data to be entered in the rate- and state- friction analysis that follows this one 
%  (data have already been already windowed, from points 1 to 4 with index specified by the DataIndex array, selected by the user using the function getpoints_velstep.m) 
%  
%  LS2: array of slip windows lengths (microns) from min_LS2 to max_LS2 with increment delta_LS2
%  steady_state: array of 1st steady-state points (microns) associated with a given moving slip window length LS2 used for the analysis 

%% How to run the code:
% To use the default input parameters:  [Time,Disp,Mu,sn,LS2,steady_state] = steadystate;

% To enter customized parameters:
% example1: change all default parameters: [Time,Disp,Mu,sn,LS2,steady_state] = steadystate(10,100,20,360);
% example2: change only deltaslip: [Time,Disp,Mu,sn,LS2,steady_state] = steadystate(10);
% example3: change deltaslip and min_LS2: [Time,Disp,Mu,sn,LS2,steady_state] = steadystate(10,100);
% example4: change deltaslip, min_LS2 and delta_LS2: [Time,Disp,Mu,sn,noise,Npointsperdisp_avs,lindetr,steady_state] = steadystate(10,100,20);

% To keep one default parameter and customize the others use [];
% example, to keep deltaslip as the default and change the others: [Time,Disp,Mu,sn,LS2,steady_state] = findsteadystate([],100,20,360)

%% Upload experimental parameters within the vel step, i.e.,DataIndex, Time, Displacement, Friction and normal stress (before step - between DataIndex 1 and 2; after step - between DataIndex 3 and 4)
close all; clc;
a = input('Have you already sliced your exp data into single velocity steps(y/n)? \n','s');
if a=='n'
    error('Slice your experiments e.g. with slicing_velsteps.m before running the program')
end
if a=='y' % import single velocity step
    [DataIndex,Time,Disp,Mu,sn]=getpoints_velstep; %call to getpoints_velsteps.m
end

%% Default checking
%Input parameters
if ((nargin < 4 || isempty(max_LS2))&&(Disp(DataIndex(4))-Disp(DataIndex(2)))<1000)
    max_LS2 = 300; %by default detrend window length after the velocity step with SL < 1000 microns = 300 microns
else
    max_LS2 = 500; %by default detrend window length after the velocity step with SL >= 1000 microns = 500 microns
end
if nargin < 3 || isempty(delta_LS2)
    delta_LS2 = 50; %by default increment in the detrend window size used to calculate the slope from one analysis to the other = 50 microns
end
if nargin < 2 || isempty(min_LS2)
    min_LS2 = 50; %by default detrend window size before the velocity step = 50 microns
end
if nargin < 1 || isempty(deltaslip)
    n=1; %counter increment in datapoints as default for the slope analysis 
         %checkslipint = Disp(DataIndex(3)+n)-Disp(DataIndex(3)); %check delta displacement right after the velocity step
else
    avgdeltadisp = abs(mean(diff(Disp(DataIndex(3):DataIndex(4))))); %average delta displacement after the velocity step
    if (deltaslip<avgdeltadisp)  %if the moving slip interval input by the user (deltaslip) is lower than 
                                 %the actual average delta displacement after the velocity step, the program by default sets
                                 %the moving slip interval to the minimum slip distance among two points and warns about the changes made
        n=1;
        warning ('input deltaslip is smaller than the equivalent slip distance between two consecutive points=');
        text = ['deltaevalslipwin has been set by default equal to the minimum delta slip = ',num2str(avgdeltadisp),' ',char(956),'m'];
        disp(text)
        %checkslipint = Disp(DataIndex(3)+n)-Disp(DataIndex(3)); %check the accuracy of the conversion in data points of deltaevalslipwin
    else
        n=round(deltaslip/avgdeltadisp); %counter increment in number of data points approximately equal to the counter step in microns
        %checkslipint = Disp(DataIndex(3)+n)-Disp(DataIndex(3)); %check the accuracy of the conversion in data points of deltaevalslipwin
    end
end 

detr = (min_LS2:delta_LS2:max_LS2)'; %intervals of slip moving windows used for evaluating steady-state

%% Calculation of slope prior to vel step between DataIndex1 and 2 (i.e., from the starting point to the one where the velocity step occurs) 

format long

tdata_bvs = Time(DataIndex(1):DataIndex(2)); %time in seconds before the velocity step
xdata_bvs = Disp(DataIndex(1):DataIndex(2)); %displacement in microns before the velocity step
ydata_bvs = Mu(DataIndex(1):DataIndex(2));  %friction coefficient before the velocity step
 
mdl = fitlm(xdata_bvs,ydata_bvs); %fit linear regression model
stud_res = mdl.Residuals.Studentized; %determination of the studentized residual to identify outliers (especially the high leverage ones) and exclude them from the slope analysis before the vel step
I_studres3 = abs(stud_res) > 3; 
outliers_studres3 = excludedata(xdata_bvs,ydata_bvs,'indices',I_studres3); %99.7% of the data are within 3 standard deviations of the mean
linreg1 = fitlm(xdata_bvs,ydata_bvs,'Exclude',outliers_studres3); %linear regression before the velocity step via removal of outliers (+-3sigma) prior to the analysis 

noise_bvs = linreg1.RMSE; %square root of the mean squared error of the linear regression (root mean square error) - proxy for the noise level of the signal before the velocity step
slope1= linreg1.Coefficients.Estimate(2); %slope before the velecity step S1
% Uncomment the line below to make the movie of the analysis
% intercept1 = linreg1.Coefficients.Estimate(1); %intersection of the linear regression

dispint = Disp(DataIndex(2))-Disp(DataIndex(1)); %displacement interval before the velocity step

deltat_bvs = mean(diff(tdata_bvs)); %delta time before the velocity step (s)
freq_bvs = 1/deltat_bvs; %average sampling frequency before the velocity step (Hz)
deltadisp_bvs = abs(mean(diff(Disp(DataIndex(1):DataIndex(2))))); %average delta displacement (microns) before the velocity step
vel_bvs = deltadisp_bvs/deltat_bvs; %average slip vel before the velocity step (microns/sec)
Npointsperdisp_bvs = freq_bvs/vel_bvs; %number of points per unit displacement before the velocity step = freq/vel = 1/deltadisp (microns^-1)

%Warnings based on linear regression analysis before the velocity step of the same synthetic velocity step 
%(i.e., same RSF paramaters a = 0.012; b1 = 0.010; b2 = 0.007; DRS_1 = 20 microns; DRS_2 = 500 microns; k'= 0.01 microns^-1)
%that differ in the superimposed noise levels, sampling frequency and slip window length of the analysis (LS1). 
%When the warning appears, it is likely that the measurement of the slope S1 might be inaccurate. 
%This might potentially affect the steady-state analysis when the calculated slope is > 5e-7 microns^-1, 
%as either the program detrends internally when not needed, as in reality S1~0, or it
%detrends of a different amount than required from the true S1 slope.

if  dispint<45
    warning('The slope measurement before the velocity step might not be accurate, hence the whole steady-state assessment: it is advised to choose either a larger slip interval for the linear regression before the vel step (LS1) or run a test at higher sampling frequency/larger step length');
    disp('flag warning 5')
    warn=['slope1 = ',num2str(slope1),' 1/',char(956),'m'];
    disp(warn)
    fprintf('\n');
end
if ((dispint<70&&dispint>=45)&&(noise_bvs<=0.0008)&&(Npointsperdisp_bvs<30))||((dispint<70&&dispint<=45)&&(noise_bvs>0.0008)&&(Npointsperdisp_bvs<3000))
    warning('The slope measurement before the velocity step might not be accurate, hence the whole steady-state assessment: it is advised to choose either a larger slip interval for the linear regression before the vel step (LS1) or run a test at higher sampling frequency/larger step length');
    disp('flag warning 4')
    warn=['slope1 = ',num2str(slope1),' 1/',char(956),'m'];
    disp(warn)
    fprintf('\n');
end 
if ((dispint<95&&dispint>=70)&&(noise_bvs<=0.000425)&&(Npointsperdisp_bvs<3))||((dispint<95&&dispint>=70)&&(noise_bvs<=0.0008&&noise_bvs>0.000425)&&(Npointsperdisp_bvs<30))||((dispint<95&&dispint>=70)&&(noise_bvs>0.0008)&&(Npointsperdisp_bvs<300))
    warning('The slope measurement before the velocity step might not be accurate, hence the whole steady-state assessment: it is advised to choose either a larger slip interval for the linear regression before the vel step (LS1) or run a test at higher sampling frequency/larger step length');
    disp('flag warning 3')
    warn=['slope1 = ',num2str(slope1),' 1/',char(956),'m'];
    disp(warn)
    fprintf('\n');
end 
if ((dispint<120&&dispint>=95)&&(noise_bvs<=0.000425&&noise_bvs>0.0003)&&(Npointsperdisp_bvs<3))||((dispint<120&&dispint>=95)&&(noise_bvs>0.000425)&&(Npointsperdisp_bvs<30))
    warning('The slope measurement before the velocity step might not be accurate, hence the whole steady-state assessment: it is advised to choose either a larger slip interval for the linear regression before the vel step (LS1) or run a test at higher sampling frequency/larger step length');
    disp('flag warning 2')
    warn=['slope1 = ',num2str(slope1),' 1/',char(956),'m'];
    disp(warn)
    fprintf('\n');
end 
if (dispint>=120&&noise_bvs>0.000425&&Npointsperdisp_bvs<3)
    warning('The slope measurement before the velocity step might not be accurate, hence the whole steady-state assessment: it is advised to choose either a larger slip interval for the linear regression before the vel step (LS1) or run a test at higher sampling frequency/larger step length');
    disp('flag warning 1')
    warn=['slope1 = ',num2str(slope1),' 1/',char(956),'m'];
    disp(warn)
    fprintf('\n');
end 

%preprocessing for data analysis: "zero-slope threshold" and removal of trends if larger than the threshold (to work with the same tolerances)

if  abs(slope1)>5e-7
    Mudetr = Mu+slope1.*(Disp(DataIndex(2))-Disp); %pre-processing the data: removing the linear trend before the velocity step using the slope prior to the velocity step S1 (only internally in the program!)
else
    Mudetr = Mu; %no need to detrend the data
end

%% Calculate the slopes after the velocity step and find the steady-state points at a given slip window length LS2
% To make a movie of the analysis, please uncomment the below lines:
% vidfile = VideoWriter('set_exp013_step3.mp4','MPEG-4');
% open(vidfile);

tdata_avs = Time(DataIndex(3):DataIndex(4)); %time in seconds after the velocity step
deltat_avs = mean(diff(tdata_avs)); %average delta time after the velocity step (s)
freq_avs = 1/deltat_avs; %average sampling frequency after the velocity step (Hz)
deltadisp_avs = abs(mean(diff(Disp(DataIndex(3):DataIndex(4))))); %delta displacement (microns) after the velocity step
vel_avs = deltadisp_avs/deltat_avs; %slip velocity after the velocity step
Npointsperdisp_avs = freq_avs/vel_avs; %number of points per unit displacement = 1/deltadisp; after the velocity step

for j = 1:length(detr) %reiterate the steady-state analysis for every detrend window length increment delta_LS2
    txt = ['Linear regression window size = ',num2str(detr(j)),' ',char(956),'m'];
    disp(txt)
    lastptevalwin = Disp(DataIndex(4))-detr(j); Nlastptevalwin = find(Disp>lastptevalwin,1)-1; Ndetr_end= length(Nlastptevalwin:DataIndex(4)); %last evaluation slip window in equivalent data points
    lastevalptdetr = Disp(DataIndex(3))+detr(j); Nlastevalptdetr = find(Disp>lastevalptdetr,1)-1; Ndetr = length(DataIndex(3):Nlastevalptdetr); %evaluation slip windows throughout the analysis in equivalent data points 
    %we made the differentiation between Ndetr_end and Ndetr since it is possible that the delta displacement might vary along the velocity step due to slight changes in velocity, 
    %Therefore, the number of data points to get the same displacement may change throughout the test.
    
    evalwin = (DataIndex(3):(DataIndex(4)- Ndetr_end))'; %overall slip window length in equivalent data points for the calculation of the slopes S2 after the velocity step  
    lastind = 1+((fix((length(evalwin)-1)/n))*n); %generalization of the last point (in data indexes) within evalwin for the slope analysis

%     Uncomment the line below for the movie of the analysis: 1st linear interpolation after velocity step with undetrended data (for figure 2)
%     linreg2 = polyfit(Disp(DataIndex(3):DataIndex(3)+Ndetr),Mu(DataIndex(3):DataIndex(3)+Ndetr),1); 

    %To quantify the mean noise after the velocity step (use ~200 microns slip window towards the end of the velocity step)
    startlinregr4noise = DataIndex(4)-round(Npointsperdisp_avs*250); endlinregr4noise = round(DataIndex(4)-Npointsperdisp_avs*50); %intervals of the linear regression model used to quantify the noise after the velocity step
    linregress4noise = fitlm(Disp(startlinregr4noise:endlinregr4noise),Mudetr(startlinregr4noise:endlinregr4noise));
    noise_avs = linregress4noise.RMSE; %square root of the mean squared error of the linear regression (root mean square error) - proxy for the noise level of the signal after the velocity step

    lindetr2 = polyfit(Disp(DataIndex(3):DataIndex(3)+Ndetr),Mudetr(DataIndex(3):DataIndex(3)+Ndetr),1); %1st linear regression after velocity step  
    slope2 = lindetr2(1); %slope S2 associated with the first regression after velocity step
    
    steadystate(j)=Disp(evalwin(1)); % initialize steady-state: the first steady-state candidate (microns) is at the first evaluation point after the velocity step (i.e., DataIndex3) 
    
    i=1; % initialize counter
       
    while (abs(slope2)>5e-7&&i<lastind) 
    %iterative loop to evaluate the first steady-state point (i.e., the first points whose slope calculated within the linear regression slip window is below or equal to the threshold 5e-7 microns^-1)
          i=i+n; %update index of n when the above condition is still met
          
          lastevalptdetr = Disp(evalwin(i))+detr(j); Nlastevalptdetr = find(Disp>lastevalptdetr,1)-1; Ndetr= length(evalwin(i):Nlastevalptdetr); 
%         uncomment the line below for the movie  
%           linreg2 = polyfit (Disp(evalwin(i):(evalwin(i)+Ndetr)),Mu(evalwin(i):(evalwin(i)+Ndetr)),1); 
          lindetr2 = polyfit (Disp(evalwin(i):(evalwin(i)+Ndetr)),Mudetr(evalwin(i):(evalwin(i)+Ndetr)),1); slope2 = lindetr2(1); 
          % update slope calculation as you move towards larger displacements along the velocity step
          steadystate(j) = Disp(evalwin(i)); %update steady-state candidate point (microns)   
          %slopetwodetr(j)=slope2(end); %to store the last slope S2 relative to the steady-state analysis
          
%           Uncomment the below lines to record a movie of the analysis
%           figure(2) % plots at different frames the interpolated slopes using the undetrended data as we move towards larger displacements
% %         along the velocity step. For computational reasons (time-efficiency) keep this as a comment while running the code with i=1.           
%           set(gcf, 'Position', get(0, 'Screensize'));
%           subplot(4,4,[1 2 5 6 9 10])
%           plot(Disp(DataIndex(1):DataIndex(4)),Mu(DataIndex(1):DataIndex(4)),'-k','Linewidth',1);hold on; %plot the friction vs. displacement series from the selected point 1 to point 4
%           plot(Disp(DataIndex(1):DataIndex(2)),slope1.*Disp(DataIndex(1):DataIndex(2))+intercept1,'-r','linewidth',1.5); hold on; %plot the linear regression line before the velocity step
%           plot(Disp(evalwin(i):(evalwin(i)+Ndetr)),linreg2(1).*Disp(evalwin(i):(evalwin(i)+Ndetr))+linreg2(2),'-r','linewidth',1.5);hold on; %plot the linear regression lines after the velocity steps as they were before the internal detrend 
%           scatter(Disp(evalwin(i)),linreg2(1).*Disp(evalwin(i))+linreg2(2),50,'g','filled','markeredgecolor','k');hold off;
%           xlabel('Displacement (\mum)'); ylabel('Friction coefficient'); 
%           ax = xticks;
%           axes=gca; axes.XAxis.Exponent = 0;  axes.YAxis.Exponent = 0;
%           F(i) = getframe(gcf); 
%           writeVideo(vidfile,F(i));        
    end
    
    if   (steadystate(j)<Disp(evalwin(lastind)))
          results = ['steadystate = ',num2str(steadystate(j)),' ',char(956),'m'];
          disp(results); %display results if steady-state condition is met before the last available point of analysis
    else
         warning('Steady-state corresponds to the last point available for the analysis and thus might not be accurate: it is advised to run another friction velocity test with longer step length')
         % it still stores a result while plotting the figure and saving the outputs, but the program gives a warning that it has been picked the last
         % point available and thus it might not be an accurate steady-state point
    end
    
    fig=figure(101); %the purpose of this plot is just to retrieve from the figure the x and y data (i.e., LS2 and first new steady-state)
    h(j)=plot(detr(j),steadystate(j));
    set(fig,'position',[-500,-500,1,1]); 
    set(fig,'visible','off'); % in order not to visualize the plot
    x(j) = get(h(j),'XData');LS2=x';
    y(j) = get(h(j),'YData');steady_state=y'; 

%     uncomment these lines to make a movie of the analysis
%     figure(2)
%     subplot(4,4,[3 4 7 8 11 12])
%     plot(steady_state,lindetr,'-ko','markersize',8,'markerfacecolor','g','markeredgecolor','k','linewidth',1);
%     xlabel ('First point at steady-state (\mum)');
%     ylabel ('Linear regression window size (\mum)');
%     yticks([detr])
%     xlim([min(ax) max(ax)])
%     ylim([detr(1)-delta_lindetr2winsize detr(end)+delta_lindetr2winsize]);
%     axes=gca; axes.XAxis.Exponent = 0;  axes.YAxis.Exponent = 0;
%     F(j) = getframe(gcf); 
%     writeVideo(vidfile,F(j));
%  
    figure(99) %plot outputs of the analysis throughout each single analysis carried out with a given linear regression window size
    plot(LS2,steady_state,'-ko','markersize',8,'markerfacecolor','k');
    ax=gca; ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0;
    xlabel ('Linear regression window length, LS_2 (\mum)');
    ylabel ('First new steady-state point (\mum)');
    xticks([detr])
    xlim([detr(1)-delta_LS2 detr(end)+delta_LS2]);
    ylim([Disp(DataIndex(3)) Disp(DataIndex(4))])
end

%   uncomment this to make a movie of the analysis
%     close(vidfile)

%display outputs in the command window
fprintf('\n');
txt1 = ['Noise level after the velocity step = ',num2str(noise_avs)];
txt2 = ['Number of points per unit displacement after the velocity step = ',num2str(Npointsperdisp_avs)]; %displays in the command window the noise level and freq/vel
disp(txt1)
disp(txt2)

%% Among the steady-state candidates for a given linear regression window length LS2, suggest the optimum pair of 1st new steady state-LS2 for linear detrend (steadybest in the program)
steadysort = sort(steady_state); %Sort results in crescent order
k=length(steady_state); %number of steady-state candidates
steadybest = []; %initialize steadybest

diffsteady = abs(diff(steady_state)); %differential of steady-state points
diffLS2 = abs(diff(LS2)); %differential of linear regression window lengths used for the analysis


%find the "steady-state" points that don't change
%significantly by varying the linear regression window length LS2, 
%and among the largest values find the first point available. If this points exists, this is the
%first new steady-state point and the associated LS2 suggested by the program

c = find(abs(diffsteady./diffLS2)<=1.5); % find optimum 1st steady state point: STEP 1)calculate the first order derivative between the differential of first steady-state points and the differential relative to LS2; 
% 1.5 is a threshold value based on experience on several velocity steps when steady-state is evaluated with deltaLS2 = 50 microns 
d = find(steady_state(c)>=steadysort(length(steadysort)-3)); %STEP 2) find the first 4 largest steady-state points satisfying the above condition

if (length(d)>0)
    steadybest = steady_state(c(d(1))); % if there is a candidate, the optimum first steady-state point is the first one on the list, i.e., the value associated with the smallest LS2 among the candidate points
    fprintf('\n'); %display combination of optimum first new steady-state point and slip window length LS2 suggested for performing the linear detrend during the subsequent RSF analysis
    txt3 = ['Optimum combination of 1st steady-state point = ',num2str(steadybest),' ',char(956),'m'];
    txt4 = ['and the associated window size to linear detrend = ',num2str(LS2(c(d(1)))),' ',char(956),'m'];
    disp(txt3);
    disp(txt4);
    fprintf('\n');
    
    figure();
    plot(LS2(c(d(1))),steady_state(c(d(1))),'-o','markersize',15,'MarkerFaceColor','r','markeredgecolor','r');hold on;
    plot(LS2(c(d(1))),steady_state(c(d(1))),'o','markersize',12,'MarkerFaceColor','k');hold on;
    plot(LS2,steady_state,'-o','markersize',8,'MarkerFaceColor','g','linewidth',1)
    axes = gca; axes.XAxis.Exponent = 0;  axes.YAxis.Exponent = 0; axes.ColorOrder = [0 0 0];
    xlabel ('Linear regression window length, LS_2 (\mum)');
    ylabel ('First new steady-state point (\mum)');
    xticks([detr])
    xlim([detr(1)-delta_LS2 detr(end)+delta_LS2]);
    ylim([Disp(DataIndex(3)) Disp(DataIndex(4))])
    close(99);
    
%     %reverse axes
%     plot(steady_state(c(d(1))),LS2(c(d(1))),'-o','markersize',15,'MarkerFaceColor','r','markeredgecolor','r');hold on;
%     plot(steady_state(c(d(1))),LS2(c(d(1))),'o','markersize',12,'MarkerFaceColor','k');hold on;
%     plot(steady_state,LS2,'-o','markersize',8,'MarkerFaceColor','g','linewidth',1)
%     axes = gca; axes.XAxis.Exponent = 0;  axes.YAxis.Exponent = 0; axes.ColorOrder = [0 0 0];
%     ylabel ('Linear regression window size (\mum)');
%     xlabel ('First steady-state point (\mum)');
%     yticks([detr])
%     ylim([detr(1)-delta_LS2 detr(end)+delta_LS2]);
%     xlim([Disp(DataIndex(3)) Disp(DataIndex(4))])
    
else
   figure();
   plot(LS2,steady_state,'-o','color','k','markersize',8,'MarkerFaceColor','g')
   axes = gca; axes.XAxis.Exponent = 0;  axes.YAxis.Exponent = 0; axes.ColorOrder = [0 0 0];
   xlabel ('Linear regression window length, LS_2 (\mum)');
   ylabel ('First new steady-state point (\mum)');
   xticks([detr])
   xlim([detr(1)-delta_LS2 detr(end)+delta_LS2]);
   ylim([Disp(DataIndex(3)) Disp(DataIndex(4))])
   close(99);
   %reverse axes
%    plot(LS2,steady_state,'-o','color','k','markersize',8,'MarkerFaceColor','g')
%    ax.ColorOrder = [0 0 0];
%    axes = gca; axes.XAxis.Exponent = 0;  axes.YAxis.Exponent = 0; axes.ColorOrder = [0 0 0];
%    ylabel ('Linear regression window size (\mum)');
%    xlabel ('First steady-state point (\mum)');
%    yticks([detr])
%    ylim([detr(1)-delta_LS2 detr(end)+delta_LS2]);
%    xlim([Disp(DataIndex(3)) Disp(DataIndex(4))])   
end    

%% Save outputs in .xlsx format
a = input('Save outputs(y/n)? \n','s');
    
%save in the workspace
Time = Time(DataIndex(1):DataIndex(4)); %slicing of data using the first and the last selected points in order 
Disp = Disp(DataIndex(1):DataIndex(4)); %to input them in the routine used to retrieve the modelled rate- and state-friction parameter values
Mu = Mu(DataIndex(1):DataIndex(4));
sn = sn(DataIndex(1):DataIndex(4));
    
if a=='n'
end
if a=='y'
    %save main outputs ("steady-state" points vs. linear detrend window) in the current folder in Excel spreadsheets
    filename=input('Save outputs: please input file name.xlsx\n','s');
    detrendtxt = ['LS2(',char(956),'m)']; steadytxt = ['1st new steady-state(',char(956),'m)'];
    hpluto = {detrendtxt   steadytxt};
    cpluto = [LS2 steady_state];
    cpluto_opt = [LS2(c(d(1))),steady_state(c(d(1)))];
    Mpluto = [hpluto;num2cell(cpluto)];
    Mpluto_opt = [hpluto;num2cell(cpluto_opt)];
    if (length(d)>0)
        xlswrite(filename,Mpluto,'all outputs')
        xlswrite(filename,Mpluto_opt,'optimum output')
    else
        xlswrite(filename,Mpluto,'all outputs')
    end 
end
end
