

clc; clear; format long;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ireadVOL = 1;

%function choices
icalkx_t    = 0; icalPSD = 0;      icalPDF = 0;
nsection    = 1; nsection_kx = 1;  ioverlap = 0; iwindow     = 0;

kloc = 1;
% information for the index range
ibe_k = 3150; iend_k = 3150; jbe_k = 1; jend_k = 1; ixwindowl= 0; ixwindowr = 0;

if ireadVOL == 0
    filepath = '/usr/local/home/glnvdb/duannas/czb58/workspace/Acoustics_Data/M6_Coldwall/TIMESERIES/';
    % Information for timeseries data file (only ntpoint is used)
    ntpoint = 1000; nskip = 0; dx_kp = 0.508545227E-04; dy_kp = 0.294736842E-04; 
    dt_sample = 5.0e-8;

    %DNS file information
    DNSIndex_k.ibe =   1181; DNSIndex_k.iend = 2350; DNSIndex_k.iskip = 1;
    DNSIndex_k.jbe =      1; DNSIndex_k.jend =  400; DNSIndex_k.jskip = 1; 
    DNSIndex_k.ibuffer = 20;

elseif ireadVOL ~= 0
    filepath = '/usr/local/home/czb58/duannas/duanl/Acoustics/M8_Sandia/run4_3200x500x600/TIMESERIES/';
    DNSIndex_k = ReadDNS_index_kplane(filepath)
    
    file_be = 228000; file_end = 228000; file_skip = 1000;
    num_file = (file_end - file_be)/file_skip + 1;
    ntpoint_total = num_file*DNSIndex_k.ntpoint; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This is the read function
%check file ReadForspectcorr.m to look at it
ReadForspectcorr




% Calculating PSD
if icalPSD == 1
    spect1dave = 0; tmean1 = 0; tmeanenergy1 = 0; tmean = 0; tmeanenergy = 0;
    numpt = 0;
    for i=1:nxpoint_k
        for j=1:nypoint_k
            spect1dtmp = fft(buffer(:,j,i,1,1));
            tmean = tmean + sum(buffer(:,j,i,1,1));
            numpt = numpt + 1;
        end
    end
    
    tmean = tmean/(ntpoint_total*numpt);
    
    var = ifft(spect1dtmp);
    
    spect1d = real(spect1dtmp)/(ntpoint_total*numpt);
    hold off
    
end % end icalPSD.eq.1

% Calculating kx-t spectrum


if icalkx_t == 1
    fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf('Calculating kx-t spectrum\n');
    
    spect2dave(ntpoint_total,nxpoint_k) = 0; tmean1 = 0; tmeanenergy1 = 0; tmean = 0; tmeanenergy = 0;
    for j=1:nypoint_k
        buffer_tmp(1:ntpoint_total,1:nxpoint_k) = buffer(:,j,:,1,1);
        size(buffer_tmp)

        spect2dtmp(:,:) = fft2(buffer_tmp)/(ntpoint_total*nxpoint_k);
       % for n=1:ntpoint_total
       %    for i=1:nxpoint_k
       %        spect2dave(n,i) = spect2dave(n,i) + spect2dtmp(n,i)*conj(spect2dtmp(n,i));
       %    end
       % end
        tmean1 = tmean1 + sum(buffer(:,j,:,1,1),1);
        tmean = tmean + sum(tmean1,3);
        tmeanenergy1 = tmeanenergy1 + sum(buffer(:,j,:,1,1).^2,1);
        tmeanenergy  = tmeanenergy + sum(tmeanenergy1,3);
    end
    %spect2dave = spect2dave/nypoint_k;
    
    
    tmean = tmean/(ntpoint_total*nxpoint_k*nypoint_k);
    tmeanenergy = tmeanenergy/(ntpoint_total*nxpoint_k*nypoint_k);
    tmeanenergy = tmeanenergy - tmean.^2;
    
 
    
    
    outfname = ['kx_t.dat'];
    outvarname = ['variables = "spectrum" \n'];
    fid = fopen(outfname,'w');
    fprintf(fid,outvarname);
    fprintf(fid,'zone i=%g, j=%g\n',ntpoint_total,nxpoint_k);
    for i=1:nxpoint_kdnum
        for n=1:ntpoint_total
            %fprintf(fid,'%17.10e',spect2dave(n,i));
            fprintf(fid,'%17.10e ',real(spect2dtmp(n,i)));
        end
    end
        
    
    outfname = ['p.dat'];
    outvarname = ['variables = "p" \n'];
    fid = fopen(outfname,'w');
    fprintf(fid,outvarname);
    fprintf(fid,'zone i=%g, j=%g, datapacking=point\n',ntpoint_total,nxpoint_k);
    for i=1:nxpoint_k
        for n=1:ntpoint_total
            fprintf(fid,'%17.10e',buffer(n,1,i,1,1));
        end
    end
    
    
end % end icalkx_t.eq.1


if icalPDF %calculate the PDF (uses Gaussian distribution)
    


%plot for each (i,j) point

for dset=1:TSkplane.dnum %here for easier modification
    
    hold on
    counter=1; %this is to help with building a legend
    for i=1:iend_k-ibe_k+1 
        for j=1:jend_k-jbe_k+1
            for k=1:1 %here for easier modification
            
 
                
data=buffer(:,j,i,k,dset)-mean(buffer(:,j,i,k,dset)); %x' =x - Mean of x

%sort data for plotting
data=sort(data);

%non dimensonalize with RMS
RMS=rms(data);
data=data/RMS;
    
%find data for PDF
MEAN=mean(data);
STD=std(data);


%calculate PDF
%PDFvalue=pdf('normal',data,MEAN,STD); Matlab built in, need Stat Toolbox
PDFvalue(:,i,j)=pdfG(MEAN,STD,data);


%Trapizoidal method to check function, should have intergral be about 1
% inter=0;
% for i=1:ntpoint-1
%     inter=inter+(PDFvalue(i)+PDFvalue(i+1))/2*(data(i+1)-data(i));
% end


%plotting function
plot(data,PDFvalue(:,i,j))
iloc=i+ibe_k-1;
jloc=j+jbe_k-1;
legendLabel{counter}=sprintf(strcat('i=', '%u',' j=', '%u  '),iloc,jloc);
counter=counter+1;
            end %k loop
        end %j loop
    end %i loop
    
%adding plot information
titlename=strcat({'PDFs of '},TSkplane.dname(dset),'''');

xlabel(strcat(TSkplane.dname(dset),'''/','$$\sqrt{',TSkplane.dname(dset),'''','^ 2}$$'),'Interpreter','latex')

ylabel('Probability')
title(titlename)
legend(legendLabel)
hold off




%Average over all points

for i=1:ntpoint
    data(i)=sum(sum(buffer(i,:,:,:,dset)));
end

data=data/((iend_k-ibe_k+1)*(jend_k-jbe_k+1));

%getting prime
data=data-mean(data);

%sort data for plotting
data=sort(data);

%non dimensonalize with RMS
RMS=rms(data);
data=data/RMS;
    
%find data for PDF
MEAN=mean(data);
STD=std(data);


%calculate PDF
AvgPDF=pdfG(MEAN,STD,data);
figure
plot(data,AvgPDF)
xlabel(strcat(TSkplane.dname(dset),'''/','$$\sqrt{',TSkplane.dname(dset),'''','^ 2}$$'),'Interpreter','latex')
ylabel('Probability')

%title generation
if ibe_k == iend_k
   ititle= sprintf(strcat('Averged PDF for i=', '%u'),ibe_k);
else
   ititle=sprintf(strcat('Averged PDF for i=', '%u',' to %u'),ibe_k,iend_k);
end

if jbe_k == jend_k
   jtitle= sprintf(strcat(' j=', '%u'),ibe_k);
else
   jtitle=sprintf(strcat(' j=', '%u',' to %u'),jbe_k,jend_k);
end
titleName=strcat(ititle,jtitle,{' of '},TSkplane.dname(dset),'''');

title(titleName)



end%data set loop





end %end icalPDF






















