% This code will read flowdata volume with dimension (kdim,idim,jdim)
% It will replace the flow field (Fname_rd) with the baseflow (Fname_baseflow)
% and write a new flow field (Fname_wt). The selected region has the
% dimension: (1:kdim,ibeg:iend,1:jdim)
%%%%%%%%%%%%%%%% Explanation of Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%File_Dim = how many dimensions that a file has
%Fname_rd = name of your input file
%Fname_wt = what to name output file
%Dnum = how many varibles to read
%Dname = the name of the variables

clc; clear; format long;

%%%%%%%%%%%%%%%% START OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
File_Dim = 3;
Fname_baseflow = '../baseflow/flowdata_00000000.h5';
Fname_rd       = '../flowdata_01310000.h5';
Fname_wt       = 'flowdata_01310000.h5';

idim_baseflow = 2170; jdim_baseflow = 320; kdim_baseflow = 295;
idim_rd       = 2170; jdim_rd       = 320; kdim_rd       = 295;
Dnum=5; Dname={'p' 'T' 'u' 'v' 'w'};

% select region for replacement
ibeg = 1800;
iend = 2170;

%%%%%%%%%%%%%%%% END OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Intialize Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fsol1=InitFlowHDF5(File_Dim);
fsol2=InitFlowHDF5(File_Dim);

fsol1.dimsm=[kdim_baseflow, idim_baseflow, jdim_baseflow];
fsol2.dimsm=[kdim_rd, idim_rd, jdim_rd];

buffer2= zeros(kdim_rd, idim_rd, jdim_rd);

iexist=0;
%%%%%%%%%%%%%%% Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=1:Dnum
  
    fsol1.dname = Dname{n};
    fsol2.dname = Dname{n};
    
    fsol1.fname=Fname_baseflow;
    fsol2.fname=Fname_rd;

    buffer1=ReadHDF5(fsol1);
    buffer2=ReadHDF5(fsol2);
    
    buffer2(:,ibeg:iend,:) = buffer1(:,ibeg:iend,:);
    
    fsol2.fname = Fname_wt;
    WriteHDF5(fsol2,buffer2,iexist);
    
    iexist=1;
    
    
    
end


%read scalar
t = ReadScalar(fsol1);
WriteScalar(fsol2,t, iexist);


        

