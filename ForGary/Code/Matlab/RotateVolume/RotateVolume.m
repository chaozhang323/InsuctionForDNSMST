clc; clear; format long;

%%%%%%%%%%%%%%%%%%Purpose%%%%%%%%%%%%%%%%%%%%%%%%%
%This program is used to modify flow volumes by flips and rotations to move
%wall locations and change the orientation while preserving the origin.
%The possible choices are: flip top to bottom, rotate clockwise and flip so 
%that a bottom wall would be on the left or right side. In the last 2 cases
%the file's x and z dimensions will be swapped
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
% FilePath: This should be the location of the grid and flowdata.
% dataFileName: This is the name of the flowdata excluding the file
% exentsion (The file will still need to be in h5 format).
% grdFileName: This is the name of the grid file excluding the file
% exentsion (The file will still need to be in h5 format).
%
%irotate: This is a selection of possible modifications
%   1: rotate clockwise 
%   2: flip top to bottom
%   3: rotate counter clockwise and flip top to bottom
%
%ioutp3d: This is a selection on wether or not to output a plt file
%%%%%%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%
% flowdata_.h5: this is the modified flow file
% grid.h5: this is the modified grid file
% flwodata_.plt: this is an optional file to be used in tecplot
%%%%%%%%%%%%%%%% START OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imax_file = 400; jmax_file = 300; kmax_file = 130;

imax_rd = 490;  jmax_rd = 20; kmax_rd = 460;  %unused

ist_rd = 1; jst_rd = 1; kst_rd = 1;  %unused

FilePath = '/usr/local/home/glnvdb/duannas/czb58/tmp/ForGary/Test_AveAcoustic/data/BotWall/REST/';
dataFileName = 'flowdata_00000000';
grdFileName ='grid';
nvar_fsol = 5; fsol_dname = {'T' 'p' 'u' 'v' 'w'};


irotate=3;
ioutp3d=1;
    
    


%%%%%%%%%%%%%%%% END OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Intialize Varibles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fsol = InitFlowHDF5();
fsolO= InitFlowHDF5();

imax=imax_file;
jmax=jmax_file;
kmax=kmax_file;

fsol.dimsm=[kmax,imax,jmax];


iexist = 0;

loc={'x' 'y' 'z'};

FlowFile=strcat(FilePath,dataFileName,'.h5');
GrdFile=strcat(FilePath,grdFileName,'.h5');

%Find index locations to help with reindexing later
for n=1:nvar_fsol
   if fsol_dname{n} == 'u'
       varu=n;
   elseif fsol_dname{n} == 'w'
       varw=n;
   end
    
end
%%%%%%%%%%%%%%%% Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch irotate 
    
    case 1
        
        if ioutp3d
           Grd=zeros(imax,kmax,jmax,1);
           Vars=zeros(imax,kmax,jmax,nvar_fsol);     
        end
        
        fsolO.dimsm=[imax,kmax,jmax];
        %read in varibles
        fsol.fname=FlowFile;
        fsolO.fname=strcat(dataFileName,'RotCW','.h5');
        for n=1:nvar_fsol
            
            fsol.dname=fsol_dname{n};
            fsolO.dname=fsol_dname{n};
            buffer=ReadHDF5_3D(fsol);
            %reindex varible locations
            buffer=permute(buffer,[2,1,3]);
            buffer=flipud(buffer);
            %buffer=rot90(buffer,3);
            
            varloc=n;
            %alter velocity to account for changed coordniate system
            if fsol.dname == 'u'
                fsolO.dname = 'w';
                buffer=-buffer;
                varloc=varw;
            elseif fsol.dname == 'w'
                fsolO.dname = 'u';
                varloc=varu;
               
            end

            WriteHDF5_3D(fsolO,buffer,iexist);
            iexist=1;

            %save data for plt output
            if ioutp3d
                Vars(:,:,:,varloc)=buffer;
            end
        end
        %get time data
        time=ReadScalar(fsol); fsolO.dimsm=[imax,kmax,jmax];
        WriteScalar(fsolO,time,iexist);
        
        iexist=0;
        fsol.fname=GrdFile;
        fsolO.fname=strcat(grdFileName,'RotCW','.h5');
        %modify grid file
        for n=1:3;
            fsol.dname=loc{n};
            fsolO.dname=loc{n};
            
            buffer=ReadHDF5_3D(fsol);
            %reindex locations
            buffer=permute(buffer,[2,1,3]);
            buffer=flipud(buffer);
           % buffer=rot90(buffer,3);
            
            varloc=n;
            %changing labels and indexing for varible directions
            if fsol.dname == 'z'
                fsolO.dname = 'x';
                varloc=1;
            elseif fsol.dname == 'x'
                fsolO.dname = 'z';
                varloc=3;
                %fixing spacing
                
                
                %fix spacing in flow direction
                for k=1:kmax
                    for j=1:jmax
                        mx=buffer(1,k,j);
                        mn=buffer(imax,k,j);
                        
                        %moving spacing and preserving origin
                        buffer(:,k,j)=-(buffer(:,k,j)-mx)+mn;
                        
                    end
                end
                
                
            end
            
            
            WriteHDF5_3D(fsolO,buffer,iexist);
            iexist=1;
            
            
            
            %save data for plt output
            if ioutp3d
                Grd(:,:,:,varloc)=buffer;
            end
        end
        
        if ioutp3d
        %setting up inputs for tecplot
        tdata=mat2tecplotSetup(nvar_fsol,fsol_dname,Grd,Vars);
        mat2tecplot(tdata,strcat(dataFileName,'RotCW','.plt'));
        end
            
            
    case 2
        
        if ioutp3d
            Vars=zeros(kmax,imax,jmax,nvar_fsol);
            Grd=zeros(kmax,imax,jmax,3);
        end
        
        
        fsolO.dimsm=fsol.dimsm;
        %read in varibles
        fsol.fname=FlowFile;
        fsolO.fname=strcat(dataFileName,'FLIP','.h5');
        for n=1:nvar_fsol
            
            fsol.dname=fsol_dname{n};
            fsolO.dname=fsol_dname{n};
            buffer=ReadHDF5_3D(fsol);
            %reindex varible locations
            buffer=flipud(buffer);
            %alter velocity to account for changed coordniate system
            if fsol.dname == 'w'
                buffer=-buffer;
            end

            WriteHDF5_3D(fsolO,buffer,iexist);
            iexist=1;
            %save data for plt output
            if ioutp3d
                Vars(:,:,:,n)=buffer;
            end
        end
        %get time data
        time=ReadScalar(fsol); 
        WriteScalar(fsolO,time,iexist);
        
        
        iexist=0;
        fsol.fname=GrdFile;
        fsolO.fname=strcat(grdFileName,'FLIP','.h5');
        %modify grid file
        for n=1:3;
            fsol.dname=loc{n};
            fsolO.dname=loc{n};
            
            buffer=ReadHDF5_3D(fsol);
            %reindex locations
            buffer=flipud(buffer);
            
            %modify z indexing
            if fsol.dname == 'z'
                for i=1:imax
                    for j=1:jmax
                        mx=buffer(1,i,j);
                        mn=buffer(kmax,i,j);
                        
                        %moving spacing and preserving origin
                        buffer(:,i,j)=-(buffer(:,i,j)-mx)+mn;
                        
                    end
                end
            end
            
            WriteHDF5_3D(fsolO,buffer,iexist);
            iexist=1;
            
            
            
            %save data for plt output
            if ioutp3d
                Grd(:,:,:,n)=buffer;
            end
                
        end
        
     
    if ioutp3d
        %setting up inputs for tecplot
        tdata=mat2tecplotSetup(nvar_fsol,fsol_dname,Grd,Vars);
        mat2tecplot(tdata,strcat(dataFileName,'FLIP','.plt'));
    end
        
        
   
    case 3
        
        if ioutp3d
           Grd=zeros(imax,kmax,jmax,1);
           Vars=zeros(imax,kmax,jmax,nvar_fsol);     
        end
        
        
        fsolO.dimsm=[imax,kmax,jmax];
        %read in varibles
        fsol.fname=FlowFile;
        fsolO.fname=strcat(dataFileName,'RotCCW','.h5');
        for n=1:nvar_fsol
            
            fsol.dname=fsol_dname{n};
            fsolO.dname=fsol_dname{n};
            buffer=ReadHDF5_3D(fsol);
            %reindex varible locations
%             buffer=permute(buffer,[2,1,3]);
%             buffer=flipud(buffer);
%             buffer=fliplr(buffer);
            buffer=rot90(buffer);
            buffer=fliplr(buffer);
            
            
            varloc=n;
            %alter velocity to account for changed coordniate system
            if fsol.dname == 'u'
                fsolO.dname = 'w';
                buffer=-buffer;
                varloc=varw;
            elseif fsol.dname == 'w'
                fsolO.dname = 'u';
                varloc=varu;
                buffer=-buffer;
               
            end

            WriteHDF5_3D(fsolO,buffer,iexist);
            iexist=1;
            
            %save data for plt output
            if ioutp3d
                Vars(:,:,:,varloc)=buffer;
            end
        end
        %get time data
        time=ReadScalar(fsol); fsolO.dimsm=[imax,kmax,jmax];
        WriteScalar(fsolO,time,iexist);
        
        iexist=0;
        fsol.fname=GrdFile;
        fsolO.fname=strcat(grdFileName,'RotCCW','.h5');
        %modify grid file
        for n=1:3;
            fsol.dname=loc{n};
            fsolO.dname=loc{n};
            
            buffer=ReadHDF5_3D(fsol);
            %reindex locations
%             buffer=permute(buffer,[2,1,3]);
%             buffer=flipud(buffer);
%             buffer=fliplr(buffer);
            
            buffer=rot90(buffer);
            buffer=fliplr(buffer);

            varloc=n;
            
            
            if fsol.dname == 'z'
                fsolO.dname = 'x';
                varloc=1;
                %fix indexing keep spacing
                for i=1:imax
                    for j=1:jmax
                        mx=buffer(i,1,j);
                        mn=buffer(i,kmax,j);
                        buffer(i,:,j)=-(buffer(i,:,j)-mx)+mn;
                    end
                end

            elseif fsol.dname == 'x'
                fsolO.dname = 'z';
                varloc=3;
                for k=1:kmax
                    for j=1:jmax
                        mx=buffer(1,k,j);
                        mn=buffer(imax,k,j);
                        buffer(:,k,j)=-(buffer(:,k,j)-mx)+mn;
                    end
                end
             
                
                
            end
            
            
            WriteHDF5_3D(fsolO,buffer,iexist);
            iexist=1;
            
            
            
            %save data for plt output
            if ioutp3d
                Grd(:,:,:,varloc)=buffer;
            end
        end
        
        if ioutp3d
        %setting up inputs for tecplot
        tdata=mat2tecplotSetup(nvar_fsol,fsol_dname,Grd,Vars);
        mat2tecplot(tdata,strcat(dataFileName,'RotCCW','.plt'));
        end
end
        