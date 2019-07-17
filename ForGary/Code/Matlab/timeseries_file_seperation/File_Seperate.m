clear all; clear;
fsol=InitFlowHDF5(4);
fsolt=InitFlowHDF5(1);
%%%%%%%%%%%%%%%% START OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FilePath = '/usr/local/home/czb58/duannas/czb58/tmp/ForGary/Mach8_Timeseries/Data/';
File_beg = 228000; File_end = 229000; File_skip= 1000 ;

imax_rd   = 1151; jmax_rd   = 500; kmax_rd   =  1; tmax_rd   = 100;
imax_file = 1151; jmax_file = 500; kmax_file = 10; tmax_file = 100;
ist_rd    = 1;    jst_rd    = 1;   kst_rd    = 5;  tst_rd    = 1;

group_name = 'kplane';
nvar_fsol = 1; fsol_dname = {'p'};

% iplane dimension
% imax_rd   = 2; jmax_rd   = 500; kmax_rd   = 600; tmax_rd   = 100;
% imax_file = 2; jmax_file = 500; kmax_file = 600; tmax_file = 100;
% group_name = 'iplane';
% should add jplane
%%%%%%%%%%%%%%%% END OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% INTIALIZE VARIBLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imax = imax_rd; jmax = jmax_rd; kmax = kmax_rd; tmax = tmax_rd;
num_files = fix((File_end-File_beg)/File_skip+1);

if group_name == 'kplane'

  fsol.dnum   = nvar_fsol;
  fsol.dimsf  = [tmax_file, jmax_file, imax_file, kmax_file];
  fsol.dimsm  = [tmax, jmax, imax, kmax];
  fsol.offset = [tst_rd-1, jst_rd-1, ist_rd-1, kst_rd-1]; 
  fsol.gname  = group_name;

  buffer=zeros(tmax, jmax, imax, kmax);

elseif group_name == 'iplane'

  fsol.dnum   = nvar_fsol;
  fsol.dimsf  = [tmax_file, jmax_file, kmax_file, imax_file];
  fsol.dimsm  = [tmax, jmax, kmax, imax];
  fsol.offset = [tst_rd-1, jst_rd-1, kst_rd-1, ist_rd-1]; 
  fsol.gname  = group_name;

  buffer=zeros(tmax, jmax, kmax, imax);    
    
end

  fsolt.dnum   = 1;
  fsolt.dimsf  = tmax_file;
  fsolt.dimsm  = tmax;
  fsolt.offset = tst_rd-1;
  fsolt.dname  = 'time';

for i=1:num_files
    file_exist=0;
    group_exist=0;
    file_num = File_beg+(i-1)*File_skip;
    fname_in = sprintf(strcat(FilePath,'timeseries_%08d.h5'),file_num);
    fname_out = sprintf(strcat('timeseries_',group_name,'_%08d.h5'),file_num);
    for n=1:nvar_fsol
        fsol.fname=fname_in;
        fsol.dname = fsol_dname{n};
        buffer = ReadHDF5(fsol);
        
        fsol.fname=fname_out;
        WriteHDF5(fsol,buffer,file_exist,group_exist);
        
        file_exist=1;
        group_exist=1;
    end
    fsolt.fname = fname_in;
    time = ReadHDF5(fsolt);
    
    fsolt.fname = fname_out;
    WriteHDF5(fsolt,time,file_exist);
    
end

