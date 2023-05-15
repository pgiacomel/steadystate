function [DataIndex,Time,Disp,Mu,sn]=getpoints_velstep 
%% Function belonging to "steadystate.m"

%  This subprogram allows the user to retrieve the DataIndex and the corresponding Time,Displacement,Friction, and Normal Stress data 
%  relative to the four following points necessary for the steady-state analysis using the program findsteadystate.m:

% 1) start point for linear detrend before the vel step;
% 2) end point for linear detrend when the step change in friction occurs;
% 3) start point for linear detrend analysis following vel step
% (suggestion: where the direct effect ends)
% 4) end point for linear detrend analysis following vel step

% N.B.: Since this is a sub-program already incorporated into findsteadystate.m 
% this program doesn't need to be run

%% Import data
close all;fclose all;clc;
format long

[FileName,PathName] = uigetfile('Select data file.txt'); % GUI file get

delimiter = '\t'; % tab delimeter (see string formatting page for more detail)
startRow = 2; % First line of numerical data in file
formatSpec = '%f%f%f%f'; % Line format spec in data array 
fileID = fopen([PathName FileName],'r'); % Open sesame
Rawdata = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);Rawdata=cell2mat(Rawdata);
Rawdata = Rawdata (1:end,:);

Time = Rawdata(:,1); %outputs required for RSFit3000: https://github.com/rmskarbek/RSFit3000
Disp = Rawdata(:,2); 
Mu = Rawdata(:,3);
sn = Rawdata(:,4); 

%% Plot velocity step and select points
exp.disp=Disp;  %displacement in microns
exp.mu=Mu; 

disp(FileName)
fprintf('\n');
disp('please select in sequential order the following points:');
fprintf('\n');
disp('1) start point of the analysis before the vel step');
disp('2) point upon the vel step, i.e., the point corresponding to the sharp jump in friction');
disp('3) start point following vel step - suggestion: where the direct effect ends (in a infinitely stiff machine)');
disp('4) last point of the analysis following vel step');
fprintf('\n');

input('press Enter to continue');

fig = gcf;
figure(gcf);
fig.Position = [400 30 880 620];
% initialState = setupFcn(gcf);     %calls the crosshair
getPoints(exp)                    %calls the function to select the points
dcm_obj = datacursormode(gcf);    %Returns the figure's data cursor mode object for customization
getcoordmouse                     % calls the function to get the mouse coordinates

disp ('selected points by mistake? press n to start over the analysis')
a = input('done(y/n)? \n\n','s');

    while a=='n'  % escape loop to cancel points collected so far
        close all;fclose all;clc;
        disp(FileName)
        fprintf('\n');
        disp('1st point: select start point of the analysis before the vel step');
        disp('2nd point: select point upon the vel step, i.e., the point corresponding to the sharp jump in friction');
        fprintf('\n');
        disp('3rd point: select start point following vel step - suggestion: where the direct effect ends (in a infinitely stiff machine)');
        disp('4th point: select last point of the analysis following vel step');
        fprintf('\n');
        
        fig = gcf;
        figure(gcf);
        fig.Position = [400 30 880 620];
%         initialState = setupFcn(gcf);     %calls the crosshair
        getPoints(exp)                    %calls the function to select the points
        dcm_obj = datacursormode(gcf);    %Returns the figure's data cursor mode object for customization
        getcoordmouse                     % calls the function to get the mouse coordinates        
        disp ('selected points by mistake? press n to start over the analysis')
        a = input('done(y/n)? \n\n','s');
    end
  if a=='y' %finished picking points
        c_info = getCursorInfo(dcm_obj);
        idx = extractfield(c_info,'DataIndex'); DataIndex = flip(idx)';  %extract DataIndex info
        DataIndex; %visualize DataIndex info in the command window
        fclose all;close all;
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

h=plot(x,y,'-','linewidth',1.5);   
ax=gca; ax.XAxis.Exponent = 0;
set(h, 'ButtonDownFcn', @lineClickedCallback); %%Select the points: Fires when mouse is pressed.
xlabel ('Displacement (\mum)');
ylabel ('Friction coefficient')
subtitle ('plot for handpicking the points for the steady-state analysis')
end

function lineClickedCallback(hLine, eventData) 
   pos = eventData.IntersectionPoint; % Get the clicked position from the event-data
   hFig = ancestor(hLine,'figure'); % Get the line's containing figure 
   cursorMode = datacursormode(hFig); % Get the figure's datacursormode object 
   cursorMode.SnaptoDataVertex = 'off';  % if on data cursors snap to nearest data 
   cursorMode.DisplayStyle = 'datatip'; % displays cursor information as a floating window and marker
   hDatatip = cursorMode.createDatatip(hLine); % Create a new data-tip at the clicked location: displays cursor information as a text box and marker
   set(hDatatip, 'Position',pos, 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor','r', 'Marker','o', 'HitTest','off');
   set(hDatatip,'OrientationMode','manual');
   set(hDatatip,'Orientation','bottomleft');
%    set(hDatatip,'UpdateFcn',@myfunction)
end

%% Customize labels for datatips
%NB: to reset the counter from one analysis to the other, type 'clear all'
%in the matlab command window

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

function Location = Get_Mouse_Location(~,~) % This function gets the location of the mouse pointer
    Location = get(gca, 'CurrentPoint');
    format shortG
    axesObjs = get(gcf);                                                                                
    dataObjs = axesObjs.Children; 
    xlimit = get(dataObjs, 'XLim'); ylimit = get(dataObjs, 'YLim'); %get the x and y limits of the plot
    mousepos=title(gca, ['[', num2str(Location(1,1)), ', ',num2str(Location(1,2)), ']  '],'Fontsize',10);
    %set(mousepos,'VerticalAlignment','bottom','HorizontalAlignment','center')
    set(mousepos,'position',[xlimit(1),ylimit(2)],'VerticalAlignment','bottom','HorizontalAlignment','left')
end

%% Creating the crosshair in the figure
% to use this uncomment the lines below and all the lines: initialState = setupFcn(gcf);
% in the section "Plot velocity step and select points"
% %NB:  The crosshair gives warnings in the Command Window when it can't
% %temporarily keep pace with the arrow pointer (e.g., out of figure's bounds), 
% %but it works fine either case, hence you can ignore the warning 
% 
% %The following functions are taken from the internal matlab input ginput: https://uk.mathworks.com/help/matlab/ref/ginput.html
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