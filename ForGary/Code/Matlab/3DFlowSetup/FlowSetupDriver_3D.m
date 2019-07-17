%%%%%%%%%%%%%%%% Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This program is used to extend the dimensions of flow and grid files in
%the j direction be copying the file along the j direction.
%
%It will create a data file containing x and z grid spacing
%
%It will find the i index that most closely corresponds to an x location
%on the upper surface of a wing

%%%%%%%%%%%%%%%% Explanation of Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%File_Dim = how many dimensions that a file has
%FlowFile_in = name of your input file
%FlowFile_out = what to name output file
%Dnum = how many varibles to read
%Dname = the name of the variables
%dimO = dimensions of your input file
%dimN = desired dimensions of your output file if not the same as input
%       data will be repeated until file is filled (only in j direction)
%GridFile_in = name of grid file to be used as input
%GridFile_out = name for output grid file
%xloc = x location to find the i index for
%iflow = wether or not a flowfile is being provided
%igrid = wether or not a gridfile is being provided


%%%%%%%%%%%%%%%% START OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
File_Dim = 3;

iflow=0;
FlowFile_in = '../REST/flowdata_00290000.h5';
FlowFile_out = 'flowdata_00290000.h5';

igrid=1;
GridFile_in='../REST/grid.h5';
GridFile_out='grid.h5';
xloc=0.00766191832;

Dnum=5; Dname={'p' 'T' 'u' 'v' 'w'};

idimO = 4949; jdimO = 1; kdimO = 295;
idimN = 4949; jdimN = 320; kdimN=295;

%%%%%%%%%%%%%%%% END OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Intialize Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fsolO=InitFlowHDF5(File_Dim);
fsolN=InitFlowHDF5(File_Dim);

fsolO.fname=FlowFile_in;
fsolN.fname=FlowFile_out;


fsolO.dimsm=[kdimO, idimO, jdimO];
fsolN.dimsm=[kdimN, idimN, jdimN];

bufferN= zeros(kdimN, idimN, jdimN);

grd={'x', 'y', 'z'};

iexist=0;
%%%%%%%%%%%%%%% Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if iflow
    
    for n=1:Dnum
        CPref=1;
        
        fsolO.dname = Dname{n};
        fsolN.dname = Dname{n};
        
        bufferO=ReadHDF5(fsolO);
        
        for nn = 1:jdimN
            
            bufferN(:,:,nn)=bufferO(:,:,CPref);
            
            if CPref == jdimO
                CPref=1;
            else
                CPref=CPref+1;
            end
        end
        
        WriteHDF5(fsolN,bufferN,iexist);
        
        iexist=1;
        
        
        
    end
    %read scalar
    t = ReadScalar(fsolO);
    WriteScalar(fsolN,t, iexist);
end



%Clone in j for grid file
iexist=0;
if igrid
    for n=3:-1:1
        
        CPref=1;
        fsolO.fname=GridFile_in;
        fsolN.fname=GridFile_out;
        
        fsolO.dname = grd{n};
        fsolN.dname = grd{n};
        
        bufferO=ReadHDF5(fsolO);
        
        for nn = 1:jdimN
            
            bufferN(:,:,nn)=bufferO(:,:,CPref);
            
            if CPref == jdimO
                CPref=1;
            else
                CPref=CPref+1;
            end
        end
        
        %pulling data about grids
        if n==1
            %finding delta x
            for i=1:idimO-1
                dx(i)=bufferO(1,i+1,1)-bufferO(1,i,1);
            end
            
            %find roughness i location
            irough=zp;
            while((bufferO(1,irough,1)-xloc) < 0)
                irough=irough+1;
            end
            
            if (abs(bufferO(1,irough,1)-xloc) > abs(bufferO(1,irough-1,1)-xloc))
                irough=irough-1;
            end
            xwall=bufferO(1,:,1);
            
        end
        
        if n==3 
            %find delta z
            for i=1:idimO-1
                dz(i)=bufferO(1,i+1,1)-bufferO(1,i,1);
            end
            %finding positive z
            
            zp=1;%defines where z becomes positive
            while(bufferO(1,zp,1) < 0)
                zp=zp+1;
            end
            %store zwall data
            zwall=bufferO(1,:,1);
        end
                
        
        WriteHDF5(fsolN,bufferN,iexist);
        
        iexist=1;
    end
    
   
   OutFile=fopen('dist.dat','wt');
   fprintf(OutFile,'%s\n','Variables = i, x_surf, dx, z_surf, dz');
   for i=1:idimO-1
       fprintf(OutFile, '%5.5u, %8.8e, %8.8e, %8.8e, %8.8e\n', ...
                        i,bufferO(1,i,1),dx(i),zwall(i),dz(i));
   end
   
   
    %plot of dx change
    plot(xwall(1:end-1),dx)
    title('Change of x grid spacing')
    xlabel('x')
    ylabel('\Delta x')
        
    
    
    
    
end

irough

