clc; clear; format long;
fsol = InitFlowHDF5();
fsol2=InitFlowHDF5();
%%%%%%%%%%%%%%%%%%Purpose%%%%%%%%%%%%%%%%%%%%%%%%%
%This program is used to compare the values of 2 HDF5 files and output the
%largest differnce of each varible and the location it happens at

%%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%
% imax_file, jmax_file, kmax_file: file dimensions (number of grid points)
%                              available in the streamwise, spanwise, and 
%                              wall-normal direction, respectively.
% imax_rd, jmax_rd, kmax_rd: number of grid points selected to read
%                            in the streamwise, the spanwise, and 
%                            the wall-normal directions, respectively. 
% ist_rd, jst_rd and kst_rd: index number of the 1st point to be read 
%                               in the steamwise, spanwise, and 
%                               wall-normal directions, respectively. 
% nvar_fsol: number of flow variable to read in the file. 
% fsol_dname: flow variable names for nvar_fsol variables, which could be
%             'u', 'v', 'w', 'p', 'T' and 'rho'. 
% FilePath: This should be the location to the files to be compared
% dataFileNames: These are the full file names of the 2 files to be
%                compared
% Each file contains one instantaneous snapshot of the flow field
% OutputFileName: This is the name assigned to the output file
%
%icompare: this defines the type of comparison being made
%           1: compares to similar files. They should have the same file
%           size and varible names
%           2: compares a volume file with a kplane. File dimensions should
%           be choosen to include the desired section of a kplane
%           3: compares a volume file with an iplane. File dimensions
%           should be choosen to include the desired section of a kplane
% 
%ioutall: Whether or not to output the differnce and error for all locations
%ivarchange: Whether or not to change varible names in Diff and Rrr outputs
%fsol_dnameNew: The new varible names for the HDF5 outputs. Their order
%               should correspond with the order in fsol_dname.
%
%n_kplane: number of k planes to be compared
%kindex[]: this should include all the indices of the kplanes
%maxFile: This is the max dimensions of the plane file (not currently used)
%k_planeImax/Jmax_rd: This is the max number of grid points to be read from the plane file 
%k_planeIst/Jst: This is the starting dimension to read from the plane file
%
%
%n_iplane: number of i planes to be compared
%iindex[]: this should include all the indices of the iplanes
%maxFile: This is the max dimensions of the plane file (not currently used)
%k_planeKmax/Jmax_rd: This is the max number of grid points to be read from the plane file 
%k_planeKst/Jst: This is the starting dimension to read from the plane file
%
%iTS: Wether or not the files contain time series data, currently only works for iplanes
%
%imfiles: Wether or not data is stored in multiple files. (currently only works with constant file sizes)
%         If yes then the constant section of file should be given in the 
%         data file place ie: timeseries_iplane
%
%nfiles: This is the number of files making up the plane data.
%
%%%%%%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%
%The maximum difference of each varible will be outputted to the choosen
%output file name in the format (this is only the case for icompare=1
%Varible,  MaxDiffernce, Error, location in zone selection(i=1 @ i_st)(i,j,k), Location in full file(i-g=1 @ 1)(i-g,j-g,k-g)
%
%Diff.h5
%This file contains the difference at each grid location
%
%Err.h5 
%This file contains the scaled error ( |(f1-f2)/f1)| ) of at each grid location
%%%%%%%%%%%%%%%% START OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imax_file = 1; jmax_file = 200; kmax_file = 125; tmax_file = 15200;

imax_rd = 1;  jmax_rd = 20; kmax_rd = 125; tmax_rd = 15200;

ist_rd = 1; jst_rd = 1; kst_rd = 1; tst_rd = 1;

FilePath = '';
dataFile1Name = '/usr/local/home/glnvdb/duannas/czb58/tmp/ForGary/TIMESERIES_iplane/Matlab_iplane/timeseries_iplane1550.h5';
dataFile2Name = '/usr/local/home/glnvdb/duannas/czb58/workspace/Acoustics_Data/M6_LargeSpan/TIMESERIES_HDF5/data/timeseries_iplane';
OutputFileName ='Diffs.dat';
nvar_fsol = 5; fsol_dname = {'T' 'p' 'u' 'v' 'w'};

iTS=1; %only avalible for iplane currently
icompare=3;


ioutall=0;
ivarchange=0;  fsol_dnameNew={'T' 'p' 'u' 'v' 'w'};

%%%%icompare=2  read kplane%%%%%%
n_kplane=1; kindex=[21]; 
k_planeImaxFile=21; k_planeJmaxFile=276; k_planetmaxFile=2;

k_planeImax_rd=21; k_planeJmax_rd=276;
k_planeIst=1;    k_planeJst=1;

%%%%icompare=3 read iplane%%%%%%%
n_iplane=1; iindex=[1550];
i_planeKmaxFile=125; i_planeJmaxFile=77; i_planetmaxFile=15200;

i_planeKmax=125; i_planeJmax=77; i_planetmax=15200;
i_planeKst=1;   i_planeJst=1;    i_planetst=1;
%%%%%when iTS=1 %%%%%
    imfiles=1; 
    nfiles=1;
    
    


%%%%%%%%%%%%%%%% END OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Intialize Varibles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imax=imax_rd+ist_rd-1;
jmax=jmax_rd+jst_rd-1;
kmax=kmax_rd+kst_rd-1;
tmax=tmax_rd+tst_rd-1;

iexist=0;
gexist=0;


MVal=0;

fsol.dnum = nvar_fsol;

dataFile1Name=strcat(FilePath,dataFile1Name);
dataFile2Name=strcat(FilePath,dataFile2Name);

% read timeseries volume file
fsol.dimsf = [kmax_file imax_file jmax_file];

if imfiles == 0
    nfiles=1;
end
    
%%%%%%%%%%%%%%%%%%% Calculating Differnces %%%%%%%%%%%%%%%%%%%%%%%%%%%

if iTS
    
    switch icompare
        
        case 1
        case 2
        case 3
            %inputing base file data
            fsol.fname = dataFile1Name;

            fsol.dimsm = [tmax_rd jmax_rd kmax_rd];
            
            File2=dataFile2Name;
            for n=1:n_iplane %for each iplane
                for counter=1:nfiles %for each file  making up iplane
               
                    %should create file name making up iplane in case of
                    %multiple files
                    if imfiles
                        dataFile2Name=sprintf(strcat(File2,'%04u_j%04u-%04u.h5'),iindex(n),jst_rd,i_planeJmax);
                    end
                
                
                    %Inputing info on k_plane file

                    % fsol2.dimsf = [i_planeKmaxFile i_planeJmaxFile
                    % i_planetmaxFile]; not used
                    fsol2.dimsm = [tmax_rd,i_planeJmax,i_planeKmax,1];
                    fsol2.offset= [i_planetst-1,i_planeJst-1,i_planeKst-1,0];
                    fsol2.rank=4;
                    fsol2.count = [1 1 1 1];
                    fsol2.stride = [1 1 1 1];
                    fsol2.offset = [0 0 0 0];
                    %setting proper offset for volume file
                    fsol.offset = [kst_rd-1 jst_rd-1 tst_rd-1];


                    %defining group name and opening output file

                    fsol2.gname='iplane';
                    OutFile=fopen(sprintf(strcat('Diff_iplane%04u_j%04u-%04u'),iindex(n),jst_rd,jmax),'wt');

                    fprintf(OutFile,'%s\n','Varible  MaxDifference   error           t     j     k     t-g   j-g   k-g');

                    for nn=1:fsol.dnum
                        fsol2.gname='iplane';
                        fsol2.fname= dataFile2Name; %


                        fsol.dname=fsol_dname{nn};
                        fsol2.dname=fsol_dname{nn};


                        buffer1=squeeze(ReadHDF5_3D(fsol));
                        buffer2=ReadHDF5_3D(fsol2);

                       %find max difference and location
                        Diff=abs(buffer1-buffer2);
                        err=Diff./buffer1;

                        Mval=Diff(1,1,1);
                        tloc=1;
                        jloc=1;
                        kloc=1;
                        for t=1:tmax_rd
                            for j=1:jmax_rd
                                for k=1:kmax_rd
                                    if Diff(t,j,k) > Mval
                                        tloc=t;
                                        jloc=j;
                                        kloc=k;
                                        Mval=Diff(t,j,k);
                                    end
                                end
                            end
                        end

                        %global locations
                        tlocg=tloc+tst_rd-1;
                        jlocg=tloc+jst_rd-1;
                        klocg=kloc+kst_rd-1;



                        if MVal==0 %if there are no differences
                            fprintf(OutFile,'%-8s,%s\n',fsol_dname{nn},'No Differences Found');%write Data
                        else

                            fprintf(OutFile,'%-8s,%8.8e ,%8.8e ,%04u ,%04u ,%04u ,%04u ,%04u ,%04u  \n'...
                                   ,fsol_dname{nn},MVal,err(tloc,jloc,kloc),tloc,jloc,kloc,tlocg ...
                                   ,jlocg,klocg);%write Data
                        end

                        if ioutall
                            if ivarchange
                               fsol2.dname=fsol_dnameNew{nn};
                            end
                            %writing output file name
                            outfilename=sprintf(strcat('iplane%04u_j%04u-%04u'),iindex(n),jst_rd,jmax);
                            fsol2.fname=strcat('Diff_',outfilename,'.h5');
                            fsol2.gname='/';
                            WriteHDF5_3D(fsol2,Diff,iexist);

                            fsol2.fname=strcat('Err_',outfilename,'.h5');
                            WriteHDF5_3D(fsol2,err,iexist);
                        end

                        iexist=1;

                    end
                    fclose(OutFile);
                    
                    iexist=0;
                    jst_rd=jst_rd+jmax_rd;
                    jmax=jmax_rd+jst_rd-1;
                    
                end
                
                
            
            end

            
            
            
            
        otherwise
            fprint('Unknown comparison choosen');
            
    end
    
    
    
    
else

if icompare==1;
    
%define file read and loops based on if full file needed
if ioutall 
    fsol.dimsm = [kmax_file imax_file jmax_file]; 
    %comparision index start
    ist=ist_rd;
    jst=jst_rd;
    kst=kst_rd;
else
    fsol.dimsm = [kmax_rd imax_rd jmax_rd];           
    fsol.offset = [kst_rd-1 ist_rd-1 jst_rd-1]; 
    %comaprision index start
    ist=1;
    jst=1;
    kst=1;
end

                   


    
OutFile=fopen(OutputFileName,'wt');
fprintf(OutFile,'%s\n','Varible  MaxDifference   error           i     j     k     i-g   j-g   k-g');
  
for n=1:fsol.dnum
    fsol.dname = fsol_dname{n}; %select varible
    
    %read in the 2 data files
    fsol.fname = dataFile1Name;
    buffer1 = ReadHDF5_3D(fsol);
    
    fsol.fname = dataFile2Name;
    buffer2 = ReadHDF5_3D(fsol);
    
    Diff=abs(buffer1-buffer2); %find differences
    err=Diff./abs(buffer1); %find err values
    
    %find max val 
    if jmax == 1
        jloc=1;
        for i=ist:imax
            for k=kst:kmax
                if Diff(k,i) > MVal
                        MVal=Diff(k,i);
                        iloc=i;
                        kloc=k;
                end
            end
        end 
    else
        for i=ist:imax
            for j=jst:jmax
                for k=kst:kmax
                    if Diff(k,i,j) > MVal
                        MVal=Diff(k,i,j);
                        iloc=i;
                        jloc=j;
                        kloc=k;
                    end
                end
            end
        end
        
    end
    
    %fixing indexes to proper value
    if ioutall && MVal~=0
        iloc=iloc-ist_rd+1;
        jloc=jloc-jst_rd+1;
        kloc=kloc-kst_rd+1;
    end
    
    ilocg=iloc+ist_rd-1;
    jlocg=jloc+jst_rd-1;
    klocg=kloc+kst_rd-1;
    

    %output max diff and location
    
    if MVal==0 %if there are no differences
        fprintf(OutFile,'%-8s,%s\n',fsol_dname{n},'No Differences Found');%write Data
    else
        
        fprintf(OutFile,'%-8s,%8.8e ,%8.8e ,%04u ,%04u ,%04u ,%04u ,%04u ,%04u  \n'...
                       ,fsol_dname{n},MVal,err(klocg,ilocg,jlocg),iloc,jloc,kloc,ilocg ...
                       ,jlocg,klocg);%write Data
    end
    
    if ioutall && MVal~=0
        
        %Writing HDF5 files
        if ivarchange
            fsol.dname=fsol_dnameNew{n};
        end
    
        fsol.fname='Diff.h5';
        WriteHDF5_3D(fsol,Diff,iexist);
        
        fsol.fname='Err.h5';
        WriteHDF5_3D(fsol,err,iexist);
        
        iexist=1;
    end

    
    MVal=0;
end
fclose(OutFile);

elseif icompare==2
    %inputing base file data
    fsol.fname = dataFile1Name;
    fsol.dimsm = [1 imax_rd jmax_rd];
    
    for n=1:n_kplane
        %Inputing info on k_plane file
        
        fsol2.dimsf = [k_planeImaxFile k_planeJmaxFile k_planetmaxFile]; %this 2 need to be not hard coded
        fsol2.dimsm = [k_planeImax_rd k_planeJmax_rd 1];
        fsol2.offset= [k_planeIst-1,k_planeJst-1,0];
        %setting proper offset for volume file
        fsol.offset = [kindex(n)-1 ist_rd-1 jst_rd-1];
        
        %defining group name and opening output file
        
        fsol2.gname=sprintf(strcat('k%4.4d'),kindex(n));
        OutFile=fopen(strcat('Diff_',fsol2.gname),'wt');
        
        fprintf(OutFile,'%s\n','Varible  MaxDifference   error           i     j     k     i-g   j-g   k-g');
        
        for nn=1:fsol.dnum
            fsol2.fname= dataFile2Name; %
            fsol2.gname=sprintf(strcat('k%4.4d'),kindex(n));
            
            fsol.dname=fsol_dname{nn};
            fsol2.dname=fsol_dname{nn};
            
            buffer1=squeeze(ReadHDF5_3D(fsol));
            buffer2=ReadHDF5_3D(fsol2);
            
           %find max difference and location
            Diff=abs(buffer1-buffer2);
            err=Diff./buffer1;
            [MVal,ilocM]=max(Diff);
            [MVal,jloc]=max(MVal);
            iloc=ilocM(jloc);
            
            
            %global locations
            ilocg=iloc+ist_rd-1;
            jlocg=jloc+jst_rd-1;
            klocg=kindex(n);
            
            
            
            
            if MVal==0 %if there are no differences
                fprintf(OutFile,'%-8s,%s\n',fsol_dname{nn},'No Differences Found');%write Data
            else
        
                fprintf(OutFile,'%-8s,%8.8e ,%8.8e ,%04u ,%04u ,%04u ,%04u ,%04u ,%04u  \n'...
                       ,fsol_dname{nn},MVal,err(iloc,jloc),iloc+k_planeIst-1,jloc+k_planeJst-1,1,ilocg ...
                       ,jlocg,klocg);%write Data
            end
        
            if ioutall
                if ivarchange
                   fsol2.dname=fsol_dnameNew{nn};
                end
                %writing output file name
                outfilename=sprintf(strcat('k%4.4d'),kindex(n));
                fsol2.fname=strcat('Diff_',outfilename,'.h5');
                fsol2.gname='/';
                WriteHDF5_3D(fsol2,Diff,iexist,gexist);
            
                fsol2.fname=strcat('Err_',outfilename,'.h5');
                WriteHDF5_3D(fsol2,err,iexist,gexist);
            end
            
            iexist=1;
            
        end
        fclose(OutFile);
        iexist=0;
    end
            
            
    
    
    
elseif icompare==3
       %inputing base file data
    fsol.fname = dataFile1Name;
  
    fsol.dimsm = [kmax_rd 1 jmax_rd];

    
    for n=1:n_iplane
        %Inputing info on k_plane file
        
        fsol2.dimsf = [i_planeKmaxFile i_planeJmaxFile i_planetmaxFile]; 
        fsol2.dimsm = [i_planeKmax k_planeJmax_rd 1];
        fsol2.offset= [i_planeKst-1,i_planeJst-1,0];
        %setting proper offset for volume file
        fsol.offset = [kst_rd-1 iindex(n)-1 jst_rd-1];
        
        %defining group name and opening output file
        
        fsol2.gname=sprintf(strcat('i%4.4d'),iindex(n));
        OutFile=fopen(strcat('Diff_',fsol2.gname),'wt');
        
        fprintf(OutFile,'%s\n','Varible  MaxDifference   error           i     j     k     i-g   j-g   k-g');
        
        for nn=1:fsol.dnum
            fsol2.fname= dataFile2Name; %
            fsol2.gname=sprintf(strcat('i%4.4d'),iindex(n));
            
            fsol.dname=fsol_dname{nn};
            fsol2.dname=fsol_dname{nn};
            
            buffer1=squeeze(ReadHDF5_3D(fsol));
            buffer2=ReadHDF5_3D(fsol2);
            
           %find max difference and location
            Diff=abs(buffer1-buffer2);
            err=Diff./buffer1;
            [MVal,klocM]=max(Diff);
            [MVal,jloc]=max(MVal);
            kloc=klocM(jloc);
            
            
            %global locations
            ilocg=iindex(n);
            jlocg=jloc+jst_rd-1;
            klocg=kloc+kst_rd-1;
            
            
            
            
            if MVal==0 %if there are no differences
                fprintf(OutFile,'%-8s,%s\n',fsol_dname{nn},'No Differences Found');%write Data
            else
        
                fprintf(OutFile,'%-8s,%8.8e ,%8.8e ,%04u ,%04u ,%04u ,%04u ,%04u ,%04u  \n'...
                       ,fsol_dname{nn},MVal,err(kloc,jloc),1,jloc+i_planeJst-1,kloc+i_planeKst-1,ilocg ...
                       ,jlocg,klocg);%write Data
            end
        
            if ioutall
                if ivarchange
                   fsol2.dname=fsol_dnameNew{nn};
                end
                %writing output file name
                outfilename=sprintf(strcat('i%4.4d'),iindex(n));
                fsol2.fname=strcat('Diff_',outfilename,'.h5');
                fsol2.gname='/';
                WriteHDF5_3D(fsol2,Diff,iexist);
            
                fsol2.fname=strcat('Err_',outfilename,'.h5');
                WriteHDF5_3D(fsol2,err,iexist);
            end
            
            iexist=1;
            
        end
        fclose(OutFile);
        iexist=0;
    end
    
end
    
end %end Timeseries option
clear n ;
