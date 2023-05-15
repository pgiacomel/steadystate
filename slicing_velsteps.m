%% PROGRAM TO SLICE THE SINGLE VELOCITY STEPS FROM RAW DATA IN .txt FORMAT  

clearvars; clc; close
%% import raw data
format long
startRow = input ("Enter the row number of the first numeric row in file .txt\n",'s');
delimiter= input('Enter the type of delimiter between the columns in file (e.g., tab= \\t, comma= ,)\n','s'); %(see string formatting page for more detail)
formatSpec = input('Input as many %f as the number of columns (e.g. 3 columns %f%f%f)\n','s'); %Line format spec in data array 

format long
disp('load raw data.txt');
[FileName,PathName] = uigetfile('Select data file.txt'); %GUI file get

fileID = fopen([PathName FileName],'r'); %Open sesame
Rawdata = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);Rawdata=cell2mat(Rawdata);
Rawdata = Rawdata (1:end,:);

disp(FileName)
fprintf('\n');

%extracting columns
a= input('Enter column number for time\n','s');
a = str2double(a); %convert the character array into a number
time=Rawdata(:,a);
a= input('Is time in seconds?(y/n)?\n','s');  
if a=='n'
    aa= input('Enter divisor factor to obtain seconds\n','s');
    aa = str2double(aa);
    time = time./aa;
else
end
b= input('Enter column number for slip\n','s');
b = str2double(b); %convert the character array into a number
slip=Rawdata(:,b);
bb = input('Zero first slip data?(y/n)?\n','s');
if bb=='y'
   slip = slip-slip(1);
else
end    
bbb= input('Is slip in microns?(y/n)?\n','s');  
if bbb=='n'
    bbbb= input('Enter multiplication coefficient to obtain microns\n','s');
    bbbb = str2double(bbbb);
    slip = slip.*bbbb;
else
end
c = input('Does the file contain the column for friction?(y/n)?\n','s');
if c=='y'
    c = input('Enter column number for friction\n','s');
    c = str2double(c);
    fric = Rawdata(:,c); 
    e= input('Working with fluid pressures?(y/n)?\n','s');
    if e=='y'
        ee = input('Enter column number for fluid pressure\n','s');
        ee = str2double(ee);
        Pf = Rawdata(:,ee);
        ee = input('Is fluid pressure in MPa?(y/n)?\n','s');
        if ee=='n'
            eee= input('Enter divisor factor to obtain MPa\n','s');
            eee = str2double(eee);
            Pf = Pf./eee;
        else
        end     
        f = input('Enter column number for normal stress\n','s');
        f = str2double(f);
        sn = Rawdata(:,f);
        f = input('Is normal stress in MPa?(y/n)?\n','s');
        if f=='n'
            ff= input('Enter divisor factor to obtain MPa\n','s');
            ff = str2double(ff);
            sn = sn./ff;
        end
      elseif e=='n' 
        f = input('Enter column for normal stress\n','s');
        f = str2double(f);
        sn = Rawdata(:,f);
        f = input('Is normal stress in MPa?(y/n)?\n','s');
        if f=='n'
            ff= input('Enter divisor factor to obtain MPa\n','s');
            ff = str2double(ff);
            sn = sn./ff;
        else
        end            
    end
elseif c=='n' 
    d= input('Enter column number for shear force\n','s');
    d = str2double(d);
    force = Rawdata(:,d);
    dd = input('Zero first force data?(y/n)?\n','s');
    if dd=='y'
        force = force-force(1);
    else
    end    
    ddd= input('Enter divisor factor to obtain shear stress (MPa) from force\n','s');
    ddd = str2double(ddd);
    tau = force./ddd;
    e= input('Working with fluid pressures?(y/n)?\n','s');
    if e=='y'
        ee = input('Enter column number for fluid pressure\n','s');
        ee = str2double(ee);
        Pf = Rawdata(:,ee);
        ee = input('Is fluid pressure in MPa?(y/n)?\n','s');
        if ee=='n'
            eee= input('Enter divisor factor to obtain MPa\n','s');
            eee = str2double(eee);
            Pf = Pf./eee;
        else
        end     
        f = input('Enter column number for normal stress\n','s');
        f = str2double(f);
        sn = Rawdata(:,f);
        f = input('Is normal stress in MPa?(y/n)?\n','s');
        if f=='n'
            ff= input('Enter divisor factor to obtain MPa\n','s');
            ff = str2double(ff);
            sn = sn./ff;
        else
        end            
        fric = tau./(sn-Pf);
    elseif e=='n' 
        f = input('Enter column for normal stress\n','s');
        f = str2double(f);
        sn = Rawdata(:,f);
        f = input('Is normal stress in MPa?(y/n)?\n','s');
        if f=='n'
            ff= input('Enter divisor factor to obtain MPa\n','s');
            ff = str2double(ff);
            sn = sn./ff;
        else
        end            
        fric = tau./sn;
    end
end

clearvars delimiter startRow formatSpec fileID ans PathName I Rawdata 

%% Routine to collect the points

%creating a data struct
exp.time=time;   %time in seconds
exp.disp=slip;   %displacement in microns
exp.mu=fric;     %friction coefficient
if e=='y'
    exp.sneff=sn-Pf; %normal stress in MPa
else
    exp.sneff=sn;    
end

clc;
disp(FileName)
fprintf('\n');
disp('for each velocity step select start and end');
file = ls('outputs.txt');       % outputs.txt contains the start and end handpicked with column 1 = time; column 2 = slip; column 3 = friction; colums 4 = normal stress
if isempty(file) == 0 
    delete outputs.txt          % delete existing output.txt file
end

input('press Enter to continue');

fig = gcf;
figure(gcf);
fig.Position = [400 30 880 620];

gca(fig); % Make sure the figure has an axis  
% initialState = setupFcn(gcf);     %calls the crosshair
getPoints(exp)                    %call the function to upload points handpicked
dcm_obj = datacursormode(gcf);    %Returns the figure's data cursor mode object for customization
getcoordmouse

disp ('selected points by mistake? press n to exit program')
a = input('done(y/n)? \n\n','s');

    if a=='n'
      clc;close all;
      disp('type clear all to reset the counter before starting a new analysis');
      return %forces the end of the program execution
    end
    if a=='y' %finished picking points      
          c_info = getCursorInfo(dcm_obj);
          idx = extractfield(c_info,'DataIndex'); DataIndex = flip(idx)';                             %extract DataIndex info
          DataIndex                                                                                   %visualize DataIndex info 
          outputs = [exp.time(DataIndex) exp.disp(DataIndex) exp.mu(DataIndex) exp.sneff(DataIndex)]; %get coordinates of the selected time(t)-slip(x)-friction(y) data point
          fprintf('[t,slip,mu,sn] = [%f, %f, %f, %f]\n', outputs')                                    %visualize coordinates in the command window
%    save('outputs.txt','outputs','-ascii','-append'); 
          disp('type clear all to reset the counter before starting a new analysis')
          save('outputs.txt','outputs','-ascii');                                                     %save file
          fclose all;
          close all;  
          fprintf('\n');
    end
    
%% save and plot sliced signal in the single vel steps 

N=numel(outputs(:,1));
disp=zeros(length(exp.disp),1);  %preallocation
mu=zeros(length(exp.disp),1);
time=zeros(length(exp.disp),1);
velstep =zeros(N/2,1);

for i=1:N/2                      %number of velocity steps
    stepstart = DataIndex(i+(i-1));
    stepend = DataIndex(2*i);
    
    time= exp.time(stepstart:stepend);  disp= exp.disp(stepstart:stepend);  mu= exp.mu(stepstart:stepend);  sneff= exp.sneff(stepstart:stepend); sneffavg = round(mean(sneff),0);
    %define the slip and friction arrays corresponding to the selected velocity step 
    text=append('velstep',num2str(i),'.txt'); %update name of different velocity steps in the different .txt files generated
    textitle =append('velstep',num2str(i));textsubtitle=append("sn'= ",num2str(sneffavg),' MPa');
    
    fileID = fopen(text,'w');
    fprintf(fileID,'%s\t%s\t%s\t%s\n','Time(s)','Displ(microns)','Friction','Effective_Normal_stress(MPa)');             
    %headers of each .txt file
    fprintf(fileID,'%f\t%f\t%f\t%f\n',[time disp mu sneff]');
    fclose(fileID);
    
    %check files
    rawdata= dlmread(text,'\t',1,0);                                    
    %rawdata is also the output of the velocity steps in the Matlab Workspace
    figure(); %check graphycally if you have cut correctly
    subplot(2,1,1) %subplot displacement vs slip
    plot(rawdata(:,2),rawdata(:,3),'-bo','markersize',0.5,'linewidth',0.1) 
    ax=gca; ax.XAxis.Exponent = 0;
    xlabel('Displacement (\mum)');ylabel('Friction');
    title(textitle,textsubtitle);
    subplot(2,1,2) %subplot displacement vs time
    plot(rawdata(:,1),rawdata(:,3),'-ko','markersize',0.5,'linewidth',0.1) 
    ax=gca; ax.XAxis.Exponent = 0;
    xlabel('Time (s)');ylabel('Friction');
end 

 %% Function to select the points for the data slicing (1 datafile = 1 velocity step)

function getPoints(varargin)
NN=numel(varargin);                                                                           %NN = number of points handpicked
for i=1:NN                             
    %t=varargin{1,i}.time;
    x=varargin{1,i}.disp;                                                                     %get the x(slip) and y(friction) coordinates for each point collected using the function outputs 
    y=varargin{1,i}.mu;
end

h=plot(x,y,'-','linewidth',1.5); 
ax=gca; ax.XAxis.Exponent = 0;
set(h, 'ButtonDownFcn',@lineClickedCallback); %Fires when mouse is pressed.
xlabel ('Displacement (\mum)');
ylabel ('Friction coefficient')
subtitle ('plot for slicing the single velocity steps from a test')
end

function lineClickedCallback(hLine, eventData)   
   pos = eventData.IntersectionPoint; % Get the clicked position from the event-data
   hFig = ancestor(hLine,'figure'); % Get the line's containing figure 
   cursorMode = datacursormode(hFig); % Get the figure's datacursormode object 
   cursorMode.SnaptoDataVertex = 'off';  % if on data cursors snap to nearest data 
   hDatatip = cursorMode.createDatatip(hLine); % Create a new data-tip at the clicked location: displays cursor information as a text box and marker  
   set(hDatatip, 'Position',pos, 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor','r', 'Marker','o', 'HitTest','off');     
   set(hDatatip,'OrientationMode','manual');
   set(hDatatip,'Orientation','bottomleft');
   set(hDatatip,'UpdateFcn',@myupdatefcn);
end

%% Customize labels for datatips  
%NB: to reset the counter from one analysis to the other, type 'clear all'
%in the matlab command window

        function [output_txt, count] = myupdatefcn(~,event_obj)
            % obj          Currently not used (empty)
            % event_obj    Handle to event object
            % output_txt   Data cursor text string (string or cell array of strings).
            persistent datatip_objects %handles to the datatips on the graph
            n = numel(datatip_objects);
            found_match = 0;
            count = 0;
            %Look for the current datatip in our list of known datatips
            while ~found_match && count < n
                count = count+1;
                found_match = (event_obj == datatip_objects{count});
            end
            if found_match
            else
                datatip_objects{n+1}=event_obj; %this datatip is new, store its handle
            end
            counter = count+1;
%            % Code from the standard 'Update Function'
            number = round(counter/2); %number of velocity step
            if mod(counter,2)==0
                   k='end'; %even number means end vel step  
            else
                   k='start';   %odd number means start vel step
            end
            output_txt = {['Step: ',num2str(number)],num2str(k)};
        end

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
    set(mousepos,'position',[xlimit(1),ylimit(2)],'VerticalAlignment','bottom','HorizontalAlignment','left')
end

%% Creating the crosshair in the figure
% to use this uncomment the lines below and all the lines: initialState = setupFcn(gcf);
% in the section "Routine to collect the points"
% %NB:  The crosshair gives warnings in the Command Window when it can't
% %temporarily keep pace with the arrow pointer (e.g., out of figure's bounds), 
% %but it works fine either case, you can ignore the warning 
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