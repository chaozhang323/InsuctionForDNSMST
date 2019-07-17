clear; close all; clc; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% built in MATLAB function to generate surface of the MATLAB logo
Z = 160*membrane(1,100); 
s=surface(Z); 
s.EdgeColor = 'none';
view(3)

I = size(Z,1); J = size(Z,2); K = 1; 
total_points = I*J*K; 

% Generate X & Y coordinates for the surface
x = linspace(1,I,I);
y = linspace(1,J,J);
[X, Y] = meshgrid(x,y); 

% Reshape to vector
X = reshape(X,total_points,1); 
Y = reshape(Y,total_points,1); 
Z = reshape(Z,total_points,1); 
vars = [X,Y,Z];

%%%%%%%%%% TecIO Setup %%%%%%%%%%%%%

tecplot_home = 'C:\Program Files\Tecplot\Tecplot 360 EX 2017 R2';
tecio_path = strcat(tecplot_home, '\bin\tecio.dll');
tecio_header_path = strcat(tecplot_home, '\include\TECIO.h');

if ~libisloaded('tecio')
    [notfound, warnings] = loadlibrary(tecio_path, tecio_header_path,...
        'alias', 'tecio');
end

libfunctionsview('tecio')

%%%%%%%%%%% Output data to .szplt data file. %%%%%%%%%%%

output_fname = 'logo.szplt'; 
dataset_title = 'MATLAB LOGO'; 
var_names = 'X, Y, Z';
file_format = 1;  % 0 - .plt;  1 - .szplt
file_type   = 0;  % 0 - grid & solution; 1 - grid only; 2 - solution only
data_type   = 2;  % 1 - single; 2 - double; ...

[isok,~,~,~,~,filehandle] = calllib('tecio', 'tecFileWriterOpen', ...
    output_fname, dataset_title, var_names, ...
    file_format, file_type, data_type, [],[] );

zname = dataset_title;
z_idx = 1;
[isok,~,~,~,~,~,~,z_idx] = calllib('tecio', 'tecZoneCreateIJK', ...
    filehandle, zname, I, J, K, [], [], [],[],0,0,0,z_idx);


for v = 1:size(vars,2)
    isok = calllib('tecio', 'tecZoneVarWriteDoubleValues', ... 
        filehandle, z_idx, v, 1, total_points, vars(:,v));
end   


calllib('tecio', 'tecFileWriterClose', filehandle); 

if libisloaded('tecio')
    unloadlibrary('tecio') 
end



