function [ tdata ] = mat2tecplotSetup(nvar,dname,Grd,Vars )
%This function sets up tdata for use in mat2tecplot, assuming data is in a
%cube form, it works for 2d data however
%
%%%%%%%%%%inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nvar = number of varibles
%dname = the a cell array containing the names of varibles
%Grd = up to a 4d matrix with the variable (x,y,z) refernce held in the last spot
%Vars = up to a 4d matrix with the variable reference in the last spot
%%%%%%%%%output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tdata = a struct used as an input for mat2tecplot
%%%%%%%%%%%%function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        loc={'x' 'y' 'z'};
        %basic inputs
        tdata.Nvar=3+nvar;
        tdata.varnames=[loc, dname];
        tdata.vformat(1:tdata.Nvar)=2;
        %setting up data
        tdata.cubes.varloc=0;
        tdata.cubes.x=Grd(:,:,:,1);
        tdata.cubes.y=Grd(:,:,:,2);
        tdata.cubes.z=Grd(:,:,:,3);
        for n=1:nvar
            tdata.cubes.v(n,:,:,:)=Vars(:,:,:,n);
        end
        

end

