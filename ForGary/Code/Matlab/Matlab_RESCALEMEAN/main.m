%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code will read the 3D (kmax,istencilsize,1) rescale_mean file
% 1. calculate the norm
% 2. output the *.plt file
% 3. convert the old format to the new format
%File_Dim = file dimensions
%FilePath_in = file path for reading
%FilePath_out = file path for writing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; format long;
File_Dim = 3;
%%%%%%%%%%%%%%%% START OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FilePath_in = '/usr/local/home/czb58/duannas/duanl/Acoustics/M8_Sandia/run3_3200x500x330/RESCALEMEAN/';
file_be = 20000; file_end = 156000; file_skip = 1000;
idim = 4; jdim = 1; kdim = 330;

iscale = 1; % uinf, rhoinf, delta0, tinf will be used if iscale = 1 (for old version data)
uinf = 1135.1; rhoinf = 0.027; delta0 = 0.035; tinf = 51.96;

%%%%%%%%%%%%%%%%%% 1.
ical_norm = 1; % whether to calculate the L2-norm

%%%%%%%%%%%%%%%%%% 2.
ioutput_plt = 1; % whether to output *.plt file, it will read the grid information for x,y,z
gridname = '/usr/local/home/czb58/duannas/duanl/Acoustics/M8_Sandia/run3_3200x500x330/REST/grid.h5';


%%%%%%%%%%%%%%%%%% 3.
iconvert_file = 0; % whether to convert the old rescalemean to the new version
FilePath_out = './';

%%%%%%%%%%%%%%%% END OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% Intialize Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dnum=6; Dname={'uave_recycle' 'wave_recycle' 'rhoave_recycle' 'tave_recycle' ...
               'nsample_rescale' 'theta_r'};
fsol_in =InitFlowHDF5(File_Dim);
fsol_out=InitFlowHDF5(File_Dim);
fsol_grid=InitFlowHDF5(File_Dim);

fsol_in.dimsm=[kdim, idim, jdim]; 
fsol_out.dimsm=[kdim, idim, jdim]; 

buffer_in= zeros(kdim, idim, jdim);


%%%%%%%%%%%%%%% Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iscale == 0
   cc = 'reading the new version data';
   disp(cc)
   uinf = 1; rhoinf = 1; delta0 = 1; 
else
   cc = 'reading the old version data';
   disp(cc)
end

if ioutput_plt == 1
   grid_dname = {'x' 'y' 'z'};
   fsol_grid.dnum = 3;
   fsol_grid.dimsm = [kdim,1,1];
   fsol_grid.gname = '/';
   buffer_grd = zeros(kdim,1,1,3);
   buffer_plt = zeros(kdim, idim, jdim);
   
   fsol_grid.fname = gridname;
   cc = strcat('reading :', fsol_grid.fname);
   disp(cc)
   for n=1:3
       fsol_grid.dname = grid_dname{n};
       buffer_grd(:,:,:,n) = ReadHDF5(fsol_grid);
   end
 
   fnum = (file_end-file_be)/file_skip + 1;
   buffer = zeros(kdim,idim,jdim,4,fnum);
    
   for n=1:fnum
        fsol_in.fname  = sprintf(strcat(FilePath_in,'rescalemean_%8.8d.h5'),file_be + (n-1)*file_skip);
        CC = strcat('reading file: ',fsol_in.fname);
        disp(CC)
       
        for nn=1:4
            fsol_in.dname = Dname{nn};
            buffer_in=ReadHDF5(fsol_in);
            if nn==1
                buffer_plt(:,:,:,1,n) = buffer_in*uinf;
            end
            if nn==2
                buffer_plt(:,:,:,2,n) = buffer_in*uinf;
            end
            if nn==3
                buffer_plt(:,:,:,3,n) = buffer_in*rhoinf;
            end
            if nn==4
                buffer_plt(:,:,:,4,n) = buffer_in*uinf*uinf;
            end        
        end
   end
   nfiles=int16((file_end-file_be)/file_skip+0.5);
   tdata=mat2tecplotSetup(4,Dname,buffer_grd,buffer_plt(:,1:1,:,:,:),nfiles);
   OutputFileName = 'rescalemean.plt';
   CC = strcat('output file: ',OutputFileName);
   disp(CC)
   mat2tecplot(tdata,OutputFileName);
   
end % end ioutput_plt == 1


if ical_norm == 1
    fnum = (file_end-file_be)/file_skip + 1;
    buffer = zeros(kdim,idim,jdim,4,fnum);
    buffer_wt = zeros(4,fnum-1);
    buffer_time = zeros(fnum);
    for n=1:fnum
        fsol_in.fname  = sprintf(strcat(FilePath_in,'rescalemean_%8.8d.h5'),file_be + (n-1)*file_skip);
        CC = strcat('reading file: ',fsol_in.fname);
        disp(CC)
       
        for nn=1:4
            fsol_in.dname = Dname{nn};
            buffer_in=ReadHDF5(fsol_in);
            if nn == 1
                buffer(:,:,:,nn,n) = buffer_in*uinf;
            elseif nn == 2
                buffer(:,:,:,nn,n) = buffer_in*uinf;
            elseif nn == 3
                buffer(:,:,:,nn,n) = buffer_in*rhoinf;
            elseif nn == 4
                buffer(:,:,:,nn,n) = buffer_in*uinf*uinf;
            end
        end
        if iscale == 0  
            fsol_in.dname = 'time';
            buffer_time(n) = ReadScalar(fsol_in);            
        end        
    end
    
    if iscale == 0
        fsol_in.fname  = sprintf(strcat(FilePath_in,'rescalemean_%8.8d.h5'),file_be); 
        fsol_in.dname = 'uinf';
        uinf = ReadScalar(fsol_in);
        fsol_in.dname = 'rhoinf';
        rhoinf = ReadScalar(fsol_in);
        fsol_in.dname = 'tinf';
        tinf = ReadScalar(fsol_in);   
        fsol_in.dname = 'delta0';
        delta0 = ReadScalar(fsol_in);    
    end
    
    filename_wt = 'rescalemean_norm.dat';
    CC = strcat('writing file: ',filename_wt);
    disp(CC);
    fwID = fopen(filename_wt,'wt');
    
    fprintf(fwID,'%s \n', 'variables=nadvance,time,tau_leto,u_norm,w_norm,rho_norm,t_norm');
    for n=2:fnum
        for nn=1:4
            buffer_tmp = buffer(:,4,1,nn,n) - buffer(:,4,1,nn,n-1);
            if nn == 1
                buffer_wt(nn,n) = norm(buffer_tmp)/uinf;  
            elseif nn == 2
                buffer_wt(nn,n) = norm(buffer_tmp)/uinf;
            elseif nn == 3
                buffer_wt(nn,n) = norm(buffer_tmp)/rhoinf;
            elseif nn == 4
                buffer_wt(nn,n) = norm(buffer_tmp)/tinf;
            end
        end
        fprintf(fwID,'%8.8d %8.8e %8.8e %8.8e %8.8e %8.8e %8.8e \n',file_be + (n-1)*file_skip,buffer_time(n),buffer_time(n)*uinf/delta0,buffer_wt(:,n) );
    end
    
 
    
end % end ical_norm == 1


if iconvert_file == 1
    fnum = (file_end-file_be)/file_skip + 1;
    for n=1:fnum
        fsol_in.fname  = sprintf(strcat(FilePath_in,'rescalemean_%8.8d.h5'),file_be + (n-1)*file_skip);
        fsol_out.fname = sprintf(strcat(FilePath_out,'rescalemean_%8.8d.h5'),file_be + (n-1)*file_skip);
        iexist=0;
    for nn=1:4
        fsol_in.dname = Dname{nn};
        fsol_out.dname = Dname{nn};
    
        buffer_in=ReadHDF5(fsol_in);
        buffer_out = zeros(kdim, idim, jdim);
        
        if nn==1
            buffer_out = buffer_in*uinf;
        end
        if nn==2
            buffer_out = buffer_in*uinf;
        end
        if nn==3
            buffer_out = buffer_in*rhoinf;
        end
        if nn==4
            buffer_out = buffer_in*uinf*uinf;
        end
        WriteHDF5(fsol_out,buffer_out,iexist);
    
        iexist=1;
    end
    fsol_in.dname = Dname{5};
    fsol_out.dname = Dname{5};
    buffer_scalar = ReadScalar(fsol_in); 
    WriteScalar(fsol_out,buffer_scalar,1);
    fsol_in.dname = Dname{6};
    fsol_out.dname = Dname{6};
    buffer_scalar = ReadScalar(fsol_in); 
    WriteScalar(fsol_out,buffer_scalar*delta0,1);
    
    
    end
    
end




        

