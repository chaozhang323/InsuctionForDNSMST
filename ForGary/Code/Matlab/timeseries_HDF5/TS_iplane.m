

clc; clear; format long;
fout=InitFlowHDF5();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filepath = '/usr/local/home/glnvdb/duannas/czb58/workspace/Acoustics_Data/M6_Coldwall/TIMESERIES/';
%filepath= '';
ireadVOL = 0;

%function choices
icalkx_t    = 0; icalPSD = 0; icalPDF=0;
%unused
nsection    = 1; nsection_kx = 1; ioverlap    = 0; iwindow     = 0;


iloc = 2156;

% information for the index range
%%%%%%%%%%%%%%%%%%%%%% for M6_Largespan  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%jbe_i =    1; jend_i = 797; kbe_i = 1; kend_i = 497; ixwindowl= 0; ixwindowr = 0;
% Information for timeseries data file (only ntpoint is used)
%ntpoint = 1000; nskip = 0; dx_kp = 0.15196998E-02; dy_kp = 0.108135169E-02; 
%dt_sample = 2.5e-7 ;
%DNS file information
%DNSIndex_i.jbe =      1; DNSIndex_i.jend =  797; DNSIndex_i.jskip = 4;
%DNSIndex_i.kbe =      1; DNSIndex_i.kend =  497; DNSIndex_i.kskip = 4; 
%DNSIndex_i.jbuffer = 20;

%%%%%%%%%%%%%%%%%%%%%% for M6_Coldwall  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jbe_i =    1; jend_i = 400; kbe_i = 1; kend_i = 460; ixwindowl= 0; ixwindowr = 0;
% Information for timeseries data file (only ntpoint is used)
ntpoint = 2; nskip = 0; dx_kp = 0.15196998E-02; dy_kp = 0.108135169E-02; 
dt_sample = 2.5e-7 ;
%DNS file information
DNSIndex_i.jbe =      1; DNSIndex_i.jend =  400; DNSIndex_i.jskip = 1;
DNSIndex_i.kbe =      1; DNSIndex_i.kend =  460; DNSIndex_i.kskip = 1; 
DNSIndex_i.jbuffer = 20;
%varible information
num_var=5; var_list={'p','T','u','v','w'};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntpoint_total = ntpoint;

if ireadVOL==0
    
    jbuffer   = DNSIndex_i.jbuffer;
    jbe_DNS   = DNSIndex_i.jbe;
    jend_DNS  = DNSIndex_i.jend;
    jskip_DNS = DNSIndex_i.jskip;
    kbe_DNS   = DNSIndex_i.kbe;
    kend_DNS  = DNSIndex_i.kend;
    kskip_DNS = DNSIndex_i.kskip;
    
    % begin, end and number of spatial points for Averaging
    if mod(jbe_i-jbe_DNS,jskip_DNS)~=0
       jbe_ave = jbe_i+jskip_DNS - mod(jbe_i-jbe_DNS,jskip_DNS);
    else
       jbe_ave = jbe_i;
    end
    nypoint_i = (jend_i-jbe_ave)/jskip_DNS + 1;
    jend_ave = jbe_ave + (nypoint_i-1)*jskip_DNS;
    
    if mod(kbe_i-kbe_DNS,kskip_DNS)~=0
       kbe_ave = kbe_k+kskip_DNS - mod(kbe_i-kbe_DNS,kskip_DNS);
    else
       kbe_ave = kbe_i;
    end
    nzpoint_i = (kend_i-kbe_ave)/kskip_DNS + 1;
    kend_ave  = kbe_ave + (nzpoint_i-1)*kskip_DNS;
    
    fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf('DNS Output timeseries spatial index range\n');
    fprintf('jbe = %4g, jend = %4g\n',jbe_DNS,jend_DNS);
    fprintf('kbe = %4g, kend = %4g\n',kbe_DNS,kend_DNS);
    fprintf('Actual Spatial Average range\n');
    fprintf('jbe_ave = %4g, jend_ave = %4g\n',jbe_ave,jend_ave);
    fprintf('kbe_ave = %4g, kend_ave = %4g\n',kbe_ave,kend_ave);
    fprintf('Number of spatial points for Averaging\n');
    fprintf('nypoint_k = %4g, nzpoint_k = %4g\n',nypoint_i,nzpoint_i);
    
    % Dimension range for TSData
    jmints_i = 1 - min(ixwindowl, (jbe_ave-jbe_DNS)/jskip_DNS);
    jmaxts_i = nypoint_i + min(ixwindowr, (jend_DNS-jend_ave)/jskip_DNS);
    javelen = jmaxts_i - jmints_i + 1;
    fprintf('streamwise Index range for kplane buffer\n');
    fprintf('jmints_i = %4g, jmaxts_i = %4g, javelen (buffer length) = %4g\n', jmints_i,jmaxts_i,javelen);
    
    jminloc = max(jbe_ave  - ixwindowl*jskip_DNS,jbe_DNS);
    jmaxloc = min(jend_ave + ixwindowr*jskip_DNS,jend_DNS);
    fprintf('Physical i-index range read TS (including xwindow)\n');
    fprintf('iminloc = %4g, imaxloc = %4g\n',jminloc,jmaxloc);
    
    % Total number of files that cover the full range of DNS output in i dir
    n_total = ceil(((jend_DNS - jbe_DNS)/jskip_DNS + 1)/jbuffer );
    
    % Determine the range of files that need to be read based on ibe_ave & iend_ave
    n_lower = 0; n_upper = 0; % File indexes that are partially read
    joffset_lower = 0; jnum_lower = 0; jnum_upper = 0;
    
    for n=1:n_total
       jbe_pfile  = jbe_DNS + ((n-1)*jbuffer)*jskip_DNS; % Could also be obtained by reading from the file name
       jend_pfile = jbe_pfile + (jbuffer-1)*jskip_DNS;  % Could also be obtained by reading from the file name
       if (jbe_pfile<=jminloc) && (jend_pfile>=jminloc)
         n_lower = n;
         joffset_lower = (jminloc - jbe_pfile)/jskip_DNS;
         jnum_lower = min( (jend_pfile - jminloc)/jskip_DNS + 1, javelen);
       end
       
       if (jbe_pfile<=jmaxloc)&&(jend_pfile>=jmaxloc)
           n_upper = n;
           if(n_upper>n_lower) 
               jnum_upper = (jmaxloc - jbe_pfile)/jskip_DNS + 1;
           end
           break
       end
    end
    
    nfile = n_upper - n_lower + 1; % Total number of files that need to be read
    
    
    nn = 0;
    for n=n_lower:n_upper
       nn = nn + 1;
       jbe_pfile  = jbe_DNS + ((n-1)*jbuffer)*jskip_DNS;
       jend_pfile = min(jbe_pfile + (jbuffer-1)*jskip_DNS, jend_DNS);
       fname_jrange(nn,:) = sprintf(strcat('j%4.4d-%4.4d'),jbe_pfile,jend_pfile);
       %sprintf(strcat('i','%4.4d-%4.4d'),ibe_pfile,iend_pfile)
    end
    
    
   
        [TSiplane(1:nfile).gname] = deal('iplane');
        [TSiplane(:).rank]  = deal(4);
    
    
    
    
    
    % Files that are fully read
    for nn=1:nfile
       TSiplane(nn).dimsf(1) = ntpoint_total;
       TSiplane(nn).dimsf(2) = jbuffer;
       TSiplane(nn).dimsf(3) = (kend_ave - kbe_ave)/kskip_DNS + 1;
       TSiplane(nn).dimsf(4) = 1;
       TSiplane(nn).offset(1) = 0;
       TSiplane(nn).offset(2) = 0;
       TSiplane(nn).offset(3) = kbe_ave/kskip_DNS-1;     %%%%%%%% editted from +1-1
       TSiplane(nn).offset(4) = 0;
    end
    fout=InitFlowHDF5();
    % Files that are partially read in j-dir
    
    %changed to 2 from 3
    TSiplane(1).dimsf(2)  = jnum_lower;
    TSiplane(1).offset(2) = joffset_lower;
    if nfile>1
       TSiplane(nfile).dimsf(2) = jnum_upper; 
    end
   
    
    for nn=1:nfile
       TSiplane(nn).dimsm(1:4)  = TSiplane(nn).dimsf(1:4);
       TSiplane(nn).block(1:4)  = TSiplane(nn).dimsf(1:4);
       TSiplane(nn).count(1:4)  = 1;
       TSiplane(nn).stride(1:4) = 1;    
    end
    
    [TSiplane.IsHSInitialized] = deal('true');
    
    
end % end ireadVOL.eq.0



 ntpoint_tmp = 0;
iexist=0;
for nnnn=1:num_var
    
   
    nnn = jmints_i;
    % reading variables
    [TSiplane(:).dname] = deal(var_list{nnnn});




    for n=1:nfile
        TSiplane(n).fname = sprintf(strcat(filepath,'timeseries_iplane%4.4d_',fname_jrange(n,:),'.h5'),iloc);  
        buffer(ntpoint_tmp+1:ntpoint_tmp+TSiplane(n).dimsf(1),nnn:(nnn+TSiplane(n).dimsf(2)-1),1:TSiplane(n).dimsf(3),1:TSiplane(n).dimsf(4))  ...
        = ReadHDF5(TSiplane(n));
        nnn = nnn + TSiplane(n).dimsf(2);
    end

    fout.fname=sprintf(strcat('timeseries_iplane%4.4d','.h5'),iloc);
    fout.dimsm=size(buffer);
    fout.dname=var_list{nnnn};
    
    WriteHDF5(fout,buffer,iexist);
    iexist=1;

end