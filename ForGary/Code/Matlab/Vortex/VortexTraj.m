clc; clear; format long;

%%%%%%%%%%%%%%%%%%Purpose%%%%%%%%%%%%%%%%%%%%%%%%%
%This is used to create a file containing the average values in the k
%direction, with the option of writing and plt file.
%%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%
% imax_file,  kmax_file: file dimensions (number of grid points)
%                              available in the streamwise, spanwise, and 
%                              wall-normal direction, respectively.
% imax_rd,  kmax_rd: number of grid points selected to read
%                            in the streamwise, the spanwise, and 
%                            the wall-normal directions, respectively. 
% ist_rd and kst_rd: index number of the 1st point to be read 
%                               in the steamwise, spanwise, and 
%                               wall-normal directions, respectively. 
% FilePath: This should be the location of the grid and flowdata.
% dataFileName: This is the name of the flowdata 
% grdFileName: This is the name of the grid file 
%
%
%ioutp3d: This is a selection on wether or not to output a plt file
%%%%%%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%
% DigitalFilter.h5: This is an output filtered file
% DigitalFilter.plt: this is an optional file to be used in tecplot
%%%%%%%%%%%%%%%% START OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%chord and beta
c=3.83095916306;
beta=739.19818;

%file name and size
fname='vortex.dat';
fsize=119;

iinter=1;%preform interpolation
fsizeN=1000; %new file size

%%%%%%%%%%%%%%%% END OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Intialize Varibles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xc=zeros(1,fsize);
ax=zeros(1,fsize);
z=zeros(1,fsize);
axN=zeros(1,fsizeN);
%%%%%%%%%%%%%%%% Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Reading File
fID=fopen(fname,'rt');

for n=1:8
    fgets(fID)
end

for n=1:fsize
    xc(n)=fscanf(fID,'%e',1);
    ax(n)=fscanf(fID,'%e',1);

end
fclose(fID);

%caclulating parameters for integration
x=xc*c;
ab=-ax/beta;

for n = 2:fsize
    z(n)=z(n-1)+trapz(x(n-1:n),ab(n-1:n));
end

fwID=fopen('Vortex_Path.dat','wt');
fprintf(fwID,'%s\n','Variables = x/c, z');

for n=1:fsize
    fprintf(fwID,'%8.8e ,%8.8e \n',xc(n),z(n));
end
fclose(fwID);


if iinter==1
    
%interpolation 
dx=(x(end)-x(1))/(fsizeN-1);
xN=x(1):dx:x(end);

% const=NatSpline([x;ax]);
% 
% loc=1;%this tracks the which spline section you are on
% for n=1:fsizeN
%     
%     if loc ~= fsize-1
%         if xN(n) > const(6,loc+1)
%             loc=loc+1;
%         end
%     end
%     
%     axN(n)=const(4,loc)*(xN(n)-const(6,loc))^3+const(3,loc)*(xN(n)-const(6,loc))^2+...
%           const(2,loc)*(xN(n)-const(6,loc))+const(1,loc);
%    
% end

axN=spline(x,ax,xN);
abN=-axN/beta;

zN=TrapInter(xN,abN);

fwID=fopen('Vortex_Path_Iterp.dat','wt');
fprintf(fwID,'%s\n','Variables = x/c, z');

for n=1:fsizeN
    fprintf(fwID,'%8.8e ,%8.8e \n',xN(n)/c,zN(n));
end
fclose(fwID);

end






    