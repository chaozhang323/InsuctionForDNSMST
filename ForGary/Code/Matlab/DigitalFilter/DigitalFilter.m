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
imax_file = 1; kmax_file = 500;

imax_rd = 1; kmax_rd = 500;  %unused

ist_rd = 1;   kst_rd = 1;  %unused

FilePath = '/usr/local/home/glnvdb/duannas/czb58/tmp/ForGary/DigitalFilter/';
dataFileName = 'spatialave_1500-1554_2D.h5';
grdFileName ='AveAcousticStat_GridMetrics.h5';




ioutp3d=1;
    
    


%%%%%%%%%%%%%%%% END OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Intialize Varibles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fsol = InitFlowHDF52D();
fsolO= InitFlowHDF5();
jmax=1;
imax=imax_file;
kmax=kmax_file;

fsol.dimsm=[kmax,imax];
fsol.gname='Stat2d';
fsolO.dimsm=[kmax,imax,1];



iexist = 0;
dataFileName=strcat(FilePath,dataFileName);
grdFileName=strcat(FilePath,grdFileName);

fsol.fname=dataFileName;
fsolO.fname='DigFilter_Stat.h5';

readVar={'rhoave' 'tave' 'uave' 'vave' 'wave' 'uv' 'uw' 'vw' 'u2' 'v2' 'w2'}; RNum=11;
writeVar={'rhoave' 'tave' 'uave' 'vave' 'wave' 'uv' 'uw' 'vw' 'uu' 'vv' 'ww'};
%%%%%%%%%%%%%%%% Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

buffer=zeros(kmax,imax,jmax,RNum);
Grd=zeros(kmax,imax,jmax,3);

        
for n=1:RNum
    fsol.dname=readVar{n};
    fsolO.dname=writeVar{n};
    
    switch n
        
        case 6
            buffer(:,:,:,n)=ReadHDF5_3D(fsol);
            buffer(:,:,:,n)=buffer(:,:,:,n)-buffer(:,:,:,3).*buffer(:,:,:,4);
            WriteHDF5_3D(fsolO,buffer(:,:,:,n),iexist);
            
            
        case 7
            buffer(:,:,:,n)=ReadHDF5_3D(fsol);
            buffer(:,:,:,n)=buffer(:,:,:,n)-buffer(:,:,:,3).*buffer(:,:,:,5);
            WriteHDF5_3D(fsolO,buffer(:,:,:,n),iexist);
            
        case 8
            buffer(:,:,:,n)=ReadHDF5_3D(fsol);
            buffer(:,:,:,n)=buffer(:,:,:,n)-buffer(:,:,:,4).*buffer(:,:,:,5);
            WriteHDF5_3D(fsolO,buffer(:,:,:,n),iexist);
            
        case 9
            buffer(:,:,:,n)=ReadHDF5_3D(fsol);
            buffer(:,:,:,n)=buffer(:,:,:,n)-buffer(:,:,:,3).*buffer(:,:,:,3);
            WriteHDF5_3D(fsolO,buffer(:,:,:,n),iexist);
            
        case 10
            buffer(:,:,:,n)=ReadHDF5_3D(fsol);
            buffer(:,:,:,n)=buffer(:,:,:,n)-buffer(:,:,:,4).*buffer(:,:,:,4);
            WriteHDF5_3D(fsolO,buffer(:,:,:,n),iexist);
            
        case 11
            buffer(:,:,:,n)=ReadHDF5_3D(fsol);
            buffer(:,:,:,n)=buffer(:,:,:,n)-buffer(:,:,:,5).*buffer(:,:,:,5);
            WriteHDF5_3D(fsolO,buffer(:,:,:,n),iexist);
            
        otherwise
            buffer(:,:,:,n)=ReadHDF5_3D(fsol);
            WriteHDF5_3D(fsolO,buffer(:,:,:,n),iexist);
            iexist=1;
    end
    
end
            
   if ioutp3d
       fsol.dname='z';
       fsol.gname='/';
       fsol.fname=grdFileName;
       
       Grd(:,:,:,3)=ReadHDF5_3D(fsol);
       
       
       tdata=mat2tecplotSetup(RNum,writeVar,Grd,buffer);
       mat2tecplot(tdata,'DigFilter_Stat.plt');
       
   end
   
   

       

            
            
            
            
            
            
            
            
            
            
            
            