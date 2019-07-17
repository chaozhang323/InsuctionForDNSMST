function [ tdata ] = mat2tecplotSetup(nvar,dname,Grd,Vars,nzone)
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

if nargin < 5
    nzone=1;
end
loc={'x' 'y' 'z'};
%basic inputs
tdata.Nvar=3+nvar;
tdata.varnames=[loc, dname];
tdata.vformat(1:tdata.Nvar)=1;
%setting up data

for nn=1:nzone
    
    tdata.cubes(nn).varloc=0;
    tdata.cubes(nn).x=Grd(:,:,:,1);
    tdata.cubes(nn).y=Grd(:,:,:,2);
    tdata.cubes(nn).z=Grd(:,:,:,3);
    for n=1:nvar
        tdata.cubes(nn).v(n,:,:,:)=Vars(:,:,:,n,nn);
    end  
end


end

