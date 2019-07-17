function [pdfvalue] = pdfG(mean,std,data)
%calculates the probility density function using a Gaussian distribution
%Takes mean, standard deviation and data set as inputs
%Outputs pdf function values at data points

n=length(data);
pdfvalue=zeros(1,n);


for i=1:n
    fx=1/(sqrt(2*pi())*std)*exp((-(data(i)-mean)^2)/(2*std^2));
    pdfvalue(i)=fx;
end

end