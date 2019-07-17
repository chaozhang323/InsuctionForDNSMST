clc; clear;

%%%%%%%%%%%%%%%% Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This program takes multiple flow data files and 1 grd file in order to
%produce a plt file with a zone for each flow file

%%%%%%%%%%%%%%%% Explanation of Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%File_Dim = how many dimensions that a file has
%FlowFileLoc = location of your flowdata files
%GridFileLoc = location of your grid file
%OutputFileName = What to name output file
%Dnum = how many varibles to read
%Dname = the name of the variables
%max = file dimensions
%st = start of data selection
%end = end of data selection
%fstart = file number you want to start reading at
%fstride = Skip these many files
%fend = stop reading files at this file number


%%%%%%%%%%%%%%%% START OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
File_Dim = 3;

FlowFileLoc='/usr/local/home/glnvdb/duannas/glnvdb/run4_792x160x295/REST/';
GridFileLoc='/usr/local/home/glnvdb/duannas/glnvdb/run4_792x160x295/REST/';
OutputFileName='Movie.plt';
Dnum=5; Dname={'p' 'T' 'u' 'v' 'w'};

imax=490;jmax=180; kmax=295;
ist = 1; jst = 10; kst=1;
iend=490; jend=10; kend=295;

fstart= 0; fstride=200000; fend=200000;

%%%%%%%%%%%%%%%% END OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Intialize Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fsol=InitFlowHDF5(File_Dim);

idimsm=iend-ist+1;
jdimsm=jend-jst+1;
kdimsm=kend-kst+1;

fsol.dimsm=[kdimsm, idimsm, jdimsm];
fsol.offset=[kst-1, ist-1, jst-1];

grdname={'x', 'y', 'z'};

nfiles=int16((fend-fstart)/fstride+0.5);
fnum=fstart;

buffer=zeros(kdimsm,idimsm,jdimsm,Dnum,nfiles);
grd=zeros(kdimsm,idimsm,jdimsm,3);

if (idimsm+ist-1 > imax || jdimsm+jst-1 > jmax || kdimsm+kst-1 > kmax)
   error='File selection outside file dimensions'
   return
end

%%%%%%%%%%%%%%% Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read grid data
fsol.fname=strcat(GridFileLoc,'grid.h5');
for n=1:3
    fsol.dname=grdname{n};
    grd(:,:,:,n)=ReadHDF5(fsol);
end

%read flow data
for n=1:nfiles
    fsol.fname=sprintf(strcat(FlowFileLoc,'flowdata_','%08u.h5'),fnum);
    
    for nn=1:Dnum
        fsol.dname=Dname{nn};
        buffer(:,:,:,nn,n)=ReadHDF5(fsol);
    end
    fnum=fnum+fstride;
end
        

tdata=mat2tecplotSetup(Dnum,Dname,grd,buffer,nfiles);

mat2tecplot(tdata,OutputFileName);


    


