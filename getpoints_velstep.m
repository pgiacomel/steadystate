function [DataIndex,Time,Disp,Mu,sn]=getpoints_velstep 
%% Subroutine of the function "steadystate.m"
% Last modified by P. Giacomel 17-Jan-2024 13:44:52 (UTC +0)

%%%  This subscript allows the user to retrieve the DataIndex and the
%  corresponding Time,Displacement,Friction, and Effective Normal Stress
%  data relative to the three following points necessary for the
%  steady-state analysis using the routine steadystate.m:

% 1) start point for linear detrend before the vel step;
% 2) end point for linear detrend when the step change in friction occurs;
% 3) end point for linear detrend analysis following vel step

%%% N.B.: Since this is a sub-routine already incorporated into steadystate.m 
% this script doesn't need to be run

%% Import data
close all;fclose all;clc;
format long

% GUI file get
[FileName,PathName] = uigetfile('Select data file.txt'); 

% tab delimeter (see string formatting page for more detail)
delimiter = '\t';
% First line of numerical data in file
startRow = 2; 
% Line format spec in data array 
formatSpec = '%f%f%f%f'; 
% Open sesame
fileID = fopen([PathName FileName],'r'); 
Rawdata = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,...
    NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);Rawdata=cell2mat(Rawdata);
Rawdata = Rawdata (1:end,:);

%outputs required for RSFit3000: https://github.com/rmskarbek/RSFit3000
Time = Rawdata(:,1); 
Disp = Rawdata(:,2); 
Mu = Rawdata(:,3);
sn = Rawdata(:,4); 

%% Plot velocity step and select points

% Displacement in microns
exp.disp=Disp;  
% Friction coefficient
exp.mu=Mu;      

disp(FileName)
fprintf('\n');
disp('please select in sequential order the following points:');
fprintf('\n');
disp('1) start point of the analysis before the velocity step');
disp('2) point upon the velocity step, i.e., the point corresponding to the sharp jump in friction');
disp('3) last point of the analysis following velocity step');
fprintf('\n');

input('press Enter to continue');

fig = gcf;
figure(gcf);
fig.Position = [400 30 880 620];

% Uncomment line below to call the crosshair
% initialState = setupFcn(gcf);     

% Calls the function to select the points
getPoints(exp)   
% Returns the figure's data cursor mode object for customization
dcm_obj = datacursormode(gcf); 
% Calls the function to get the mouse coordinates
getcoordmouse                     

disp ('selected points by mistake? press n to start over the analysis')
a = input('done(y/n)? \n\n','s');

    % Escape loop to cancel points collected so far
    while a=='n'  
        close all;fclose all;clc;
        disp(FileName)
        fprintf('\n');
        disp('1st point: select start point of the analysis before the velocity step');
        disp('2nd point: select point upon the velocity step, i.e., the point corresponding to the sharp jump in friction');
        disp('3rd point: select last point of the analysis following velocity step');
        fprintf('\n');
        
        fig = gcf;
        figure(gcf);
        fig.Position = [400 30 880 620];
        
        % Uncomment line below to call the crosshair
%         initialState = setupFcn(gcf);   
        
        %Calls the function to select the points
        getPoints(exp) 
        %Returns the figure's data cursor mode object for customization
        dcm_obj = datacursormode(gcf); 
        %Calls the function to get the mouse coordinates  
        getcoordmouse                           
        disp ('selected points by mistake? press n to start over the analysis')
        a = input('done(y/n)? \n\n','s');
    end
  %Press y when all the 3 points have been picked   
  if a=='y' 
        c_info = getCursorInfo(dcm_obj);
        % Extract DataIndex info
        idx = extractfield(c_info,'DataIndex'); DataIndex = flip(idx)'; 
        % Visualize DataIndex info in the command window
        DataIndex; 
        fclose all;close all;
        % Uncomment line below when customizing labels for datatips
%         disp ('type clear all to reset the counter before starting a new analysis')
  end
end

%% Get coordinates of click

function getPoints(varargin)
NN=numel(varargin);                                                                                           
for i=1:NN                             
    x=varargin{1,i}.disp;                                                                                        
    y=varargin{1,i}.mu;
end

% Plot friction vs. displacement data relative to the selected velocity step
h=plot(x,y,'-','linewidth',1.5); 
ax=gca; ax.XAxis.Exponent = 0;
% Select the points: Fires when mouse is pressed.
set(h, 'ButtonDownFcn', @lineClickedCallback); 
xlabel ('Displacement (\mum)');
ylabel ('Friction coefficient')
subtitle ('plot for handpicking the points for the steady-state analysis')
end

function lineClickedCallback(hLine, eventData) 
   % Get the clicked position from the event-data
   pos = eventData.IntersectionPoint;
   % Get the line's containing figure
   hFig = ancestor(hLine,'figure'); 
   % Get the figure's datacursormode object
   cursorMode = datacursormode(hFig);  
   % If on data cursors snap to nearest data 
   cursorMode.SnaptoDataVertex = 'off';  
   % Displays cursor information as a floating window and marker
   cursorMode.DisplayStyle = 'datatip'; 
   % Create a new data tip at the clicked location: displays cursor
   % information as a text box and marker
   hDatatip = cursorMode.createDatatip(hLine); 
   set(hDatatip, 'Position',pos, 'MarkerSize',8, 'MarkerEdgeColor','k',...
       'MarkerFaceColor','r', 'Marker','o', 'HitTest','off');
   set(hDatatip,'OrientationMode','manual');
   set(hDatatip,'Orientation','bottomleft');
   set(hDatatip,'FontSize',8);
   %Uncomment line below when customizing labels for datatips
%    set(hDatatip,'UpdateFcn',@myfunction) 
end

%% Customize labels for datatips
%%% Uncomment function below to customize labels for datatips
%%% Note: to reset the counter from one analysis to the other, type 'clear all'
%%% in the matlab command window

% function [output_txt, count] = myfunction(~,event_obj)
%      % obj          Currently not used (empty)
%      % event_obj    Handle to event object
%      % output_txt   Data cursor text string (string or cell array of strings).
%      persistent datatip_objects %handles to the datatips on the graph
%      n = numel(datatip_objects);
%      found_match = 0;
%      count = 0;
%      %Look for the current datatip in our list of known datatips
%      while ~found_match && count < n
%            count = count+1;
%            found_match = (event_obj == datatip_objects{count});
%      end
%      if found_match
%      else
%            datatip_objects{n+1}=event_obj; %this datatip is new, store its handle
%      end
%      counter = count+1;
% %    % Code from the standard 'Update Function'
%      output_txt = {counter};
% end

%% Get coordinates of mouse among x-y data  

function getcoordmouse        
    set(gcf, 'WindowButtonMotionFcn', @Get_Mouse_Location);
end

% This function gets the location of the mouse pointer
function Location = Get_Mouse_Location(~,~) 
    Location = get(gca, 'CurrentPoint');
    format shortG
    axesObjs = get(gcf);                                                                                
    dataObjs = axesObjs.Children; 
    %get the x and y limits of the plot
    xlimit = get(dataObjs, 'XLim'); ylimit = get(dataObjs, 'YLim'); 
    mousepos=title(gca, ['[', num2str(Location(1,1)), ', ',...
        num2str(Location(1,2)), ']  '],'Fontsize',10);
    %set(mousepos,'VerticalAlignment','bottom','HorizontalAlignment','center')
    set(mousepos,'position',[xlimit(1),ylimit(2)],'VerticalAlignment','bottom',...
        'HorizontalAlignment','left')
end

%% Creating the crosshair in the figure
% To unlock the crosshair, uncomment the functions below and all the lines:
% initialState = setupFcn(gcf);
% in the section "Plot velocity step and select points" (lines 64,91).

%%%Note:  The crosshair gives warnings in the Command Window when it can't
% temporarily keep pace with the arrow pointer (e.g., out of figure's
% bounds), %but it works fine either case, hence you can ignore the
% warning.
% 
% The following functions are taken from the internal matlab input ginput:
% https://uk.mathworks.com/help/matlab/ref/ginput.html


% function updateCrossHair(fig, crossHair)
% % update cross hair for figure.
% gap = 3; % 3 pixel view port between the crosshairs
% cp = hgconvertunits(fig, [fig.CurrentPoint 0 0], fig.Units, 'pixels', fig);
% cp = cp(1:2);
% figPos = hgconvertunits(fig, fig.Position, fig.Units, 'pixels', fig.Parent);
% figWidth = figPos(3);
% figHeight = figPos(4);
% 
% % Early return if point is outside the figure
% set(crossHair, 'Visible', 'on');
% thickness = 1; % 1 Pixel thin lines. 
% set(crossHair(1), 'Position', [0 cp(2) cp(1)-gap thickness]);
% set(crossHair(2), 'Position', [cp(1)+gap cp(2) figWidth-cp(1)-gap thickness]);
% set(crossHair(3), 'Position', [cp(1) 0 thickness cp(2)-gap]);
% set(crossHair(4), 'Position', [cp(1) cp(2)+gap thickness figHeight-cp(2)-gap]);
% end
% 
% function crossHair = createCrossHair(fig)
% % Create thin uicontrols with black backgrounds to simulate fullcrosshair pointer.
% % 1: horizontal left, 2: horizontal right, 3: vertical bottom, 4: vertical top
% 
% if isWebFigureType(fig, 'UIFigure')
%     for k = 1:4
%         crossHair(k) = uilabel(fig, 'Visible', 'off', 'BackgroundColor', [0 0 0], 'HandleVisibility', 'off'); %#ok<AGROW>
%     end
% else
%     for k = 1:4
%         crossHair(k) = uicontrol(fig, 'Style', 'text', 'Visible', 'off', 'Units', 'pixels', 'BackgroundColor', [0 0 0], 'HandleVisibility', 'off', 'HitTest', 'off'); %#ok<AGROW>
%     end
% end
% 
% end
% 
% function initialState = setupFcn(fig)
% % Create uicontrols to simulate fullcrosshair pointer.
% initialState.CrossHair = createCrossHair(fig);
% 
% % Adding this to enable automatic updating of currentpoint on the figure 
% % This function is also used to update the display of the fullcrosshair
% % pointer and make them track the currentpoint.
% set(fig,'WindowButtonMotionFcn',@(o,e) dummy()); % Add dummy so that the CurrentPoint is constantly updated
% initialState.MouseListener = addlistener(fig,'WindowMouseMotion', @(o,e) updateCrossHair(o,initialState.CrossHair));
% end
% 
% function dummy(~,~) 
% end