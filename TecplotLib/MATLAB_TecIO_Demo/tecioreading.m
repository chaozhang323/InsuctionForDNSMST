clear all; close all; clc; 

%%%%%%%%%% TecIO Setup %%%%%%%%%%%%%
tecplot_home = 'C:\Program Files\Tecplot\Tecplot 360 EX 2017 R2';
tecio_path = strcat(tecplot_home, '\bin\tecio.dll');
tecio_header_path = strcat(tecplot_home, '\include\TECIO.h');
if ~libisloaded('tecio')
    [notfound,warnings]=loadlibrary(tecio_path, tecio_header_path, ...
        'alias', 'tecio');
end

%libfunctionsview('tecio')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pass any .szplt file here
testinfile = 'logo.szplt'; 
% Pass NULL by [] to get the filehandle
[isok, ~, handle] = calllib('tecio', 'tecFileReaderOpen', testinfile, []); 

% Get number of Variables and Zones
numvars = 0;  % Defining variable
numzones = 0; % Defining variable
[isok, ~, numvars] = calllib('tecio', 'tecDataSetGetNumVars', handle,  numvars)
[isok, ~, numzones] = calllib('tecio', 'tecDataSetGetNumZones', handle,  numzones)

for zone = 1:numzones
    % Get information on what the Zone Type is. 
    type= 0; 
    [isok, ~, type] = calllib('tecio', 'tecZoneGetType', handle, zone, type)
    
    % 0 - Ordered, 1- FE line, ... 5- FE Brick...
        
    if type == 0 
        I = 0; J = 0; K = 0; 
        [isok, ~, I, J, K] = calllib('tecio', 'tecZoneGetIJK', handle, zone, I, J, K);
    end
      
    for var = 1:numvars
        numvals = 0; 
        [isok, ~, numvals] = calllib('tecio', 'tecZoneVarGetNumValues',...
            handle, zone,var, numvals);
        values = zeros(1,numvals); % Set up space for return value
        [isok, ~, values] = calllib('tecio', 'tecZoneVarGetFloatValues', ...
            handle, zone,var, 1,numvals,values);
    end
    
end

% To get a string out must create an empty cell sructure 
title = libpointer('stringPtrPtr', cell(1,1));
[isok, ~, title] = calllib('tecio', 'tecDataSetGetTitle', handle,  title)

calllib('tecio', 'tecFileReaderClose', handle) ;

if libisloaded('tecio')
    unloadlibrary('tecio') 
end
