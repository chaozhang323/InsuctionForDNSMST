clc; clear; format long;

%%%%%%%%%%%%%%%%%%Purpose%%%%%%%%%%%%%%%%%%%%%%%%%
%This creates a grid spacing based on the desired ratio and number of
%points
%%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%
% N = number of grid points
% ratio = grid stretching resolution (zL/z2) 
%%%%%%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%
% dist_geo= Geometric Distrubtion 
%%%%%%%%%%%%%%%% START OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=5;
ratio=8;
%%%%%%%%%%%%%%%% END OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Intialize Varibles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N < 1
   error='N to small for geometric dist'
   return
end


fnc= @(x) (x^(N-1)-1)/(x-1)-ratio; %function to drive to 0

%Initilize for secant method
maxiter = 1000000;
tol=1e-5;
ds1=1/ratio;
a1=1.001; a2=1.002;
f1=fnc(a1);
f2=fnc(a2);
n=0;

dist_geo=zeros(1,N);

%%%%%%%%%%%%%%%% Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%using matlab built in
alpha=fzero(fnc,a1);


%copy from the code
% while (n < maxiter && abs(f2) > tol)
%     fdf= f2*(a2-a1)/(f2-f1);
%     if abs(fdf) > 0.1
%         fdf=fdf/abs(fdf)*0.1;
%     end
%     a1=a2;
%     f1=f2;
%     a2=a2-fdf;
%     f2=fnc(a2);
%     n=n+1;
% end
% alpha=a2;

%cacluating distribution
for i=1:N
    if i==1
        dist_geo(i)=0;
    else
        dist_geo(i)=ds1*(alpha^(i-1)-1)/(alpha-1);
    end
end



