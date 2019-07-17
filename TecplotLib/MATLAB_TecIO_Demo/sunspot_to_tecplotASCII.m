clear; close all; clc; 
load sunspot.dat;

plot(sunspot(:,1), sunspot(:,2)) 

datfile = fopen('sunspot_tec.dat', 'w+');

%Save variable list information to ASCII file
fprintf(datfile, 'VARIABLES = "Year" "Number of Sunspots"\n');

%----------ASCII----------
%Zone names and line size. Each file has own zone.
% Creates a 1D line type (See Data Format Guide)
I = size(sunspot,1);
J = 1;
K = 1;
fprintf(datfile, 'ZONE T = "Sunspot" \n');
fprintf(datfile, '    I = %i, J = %i, K = %i \n',I,J,K);

% Loop through all points along a line
for u=1:I
    % Loop through all variables 
    for v=1:size(sunspot,2) 
        fprintf(datfile,'    %f ',sunspot(u,v));
    end
    fprintf(datfile,'\n');
end

fclose(datfile);

