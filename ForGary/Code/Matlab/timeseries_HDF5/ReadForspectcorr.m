if ireadVOL==0
    ntpoint_total = ntpoint;
    ibuffer   = DNSIndex_k.ibuffer;
    ibe_DNS   = DNSIndex_k.ibe;
    iend_DNS  = DNSIndex_k.iend;
    iskip_DNS = DNSIndex_k.iskip;
    jbe_DNS   = DNSIndex_k.jbe;
    jend_DNS  = DNSIndex_k.jend;
    jskip_DNS = DNSIndex_k.jskip;
    
    % begin, end and number of spatial points for Averaging
    if mod(ibe_k-ibe_DNS,iskip_DNS)~=0
       ibe_ave = ibe_k+iskip_DNS - mod(ibe_k-ibe_DNS,iskip_DNS);
    else
       ibe_ave = ibe_k;
    end
    nxpoint_k = (iend_k-ibe_ave)/iskip_DNS + 1;
    iend_ave = ibe_ave + (nxpoint_k-1)*iskip_DNS;
    
    if mod(jbe_k-jbe_DNS,jskip_DNS)~=0
       jbe_ave = jbe_k+jskip_DNS - mod(jbe_k-jbe_DNS,jskip_DNS);
    else
       jbe_ave = jbe_k;
    end
    nypoint_k = (jend_k-jbe_ave)/jskip_DNS + 1;
    jend_ave  = jbe_ave + (nypoint_k-1)*jskip_DNS;
    fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf('DNS Output timeseries spatial index range\n');
    fprintf('ibe = %4g, iend = %4g\n',ibe_DNS,iend_DNS);
    fprintf('jbe = %4g, jend = %4g\n',jbe_DNS,jend_DNS);
    fprintf('Actual Spatial Average range\n');
    fprintf('ibe_ave = %4g, iend_ave = %4g\n',ibe_ave,iend_ave);
    fprintf('jbe_ave = %4g, jend_ave = %4g\n',jbe_ave,jend_ave);
    fprintf('Number of spatial points for Averaging\n');
    fprintf('nxpoint_k = %4g, nypoint_k = %4g\n',nxpoint_k,nypoint_k);
    
    % Dimension range for TSData
    imints_k = 1 - min(ixwindowl, (ibe_ave-ibe_DNS)/iskip_DNS);
    imaxts_k = nxpoint_k + min(ixwindowr, (iend_DNS-iend_ave)/iskip_DNS);
    iavelen = imaxts_k - imints_k + 1;
    fprintf('streamwise Index range for kplane buffer\n');
    fprintf('imints_k = %4g, imaxts_k = %4g, iavelen (buffer length) = %4g\n', imints_k,imaxts_k,iavelen);
    
    iminloc = max(ibe_ave  - ixwindowl*iskip_DNS,ibe_DNS);
    imaxloc = min(iend_ave + ixwindowr*iskip_DNS,iend_DNS);
    fprintf('Physical i-index range read TS (including xwindow)\n');
    fprintf('iminloc = %4g, imaxloc = %4g\n',iminloc,imaxloc);
    
    % Total number of files that cover the full range of DNS output in i dir
    n_total = ceil(((iend_DNS - ibe_DNS)/iskip_DNS + 1)/ibuffer );
    
    % Determine the range of files that need to be read based on ibe_ave & iend_ave
    n_lower = 0; n_upper = 0; % File indexes that are partially read
    ioffset_lower = 0; inum_lower = 0; inum_upper = 0;
    
    for n=1:n_total
       ibe_pfile  = ibe_DNS + ((n-1)*ibuffer)*iskip_DNS; % Could also be obtained by reading from the file name
       iend_pfile = ibe_pfile + (ibuffer-1)*iskip_DNS;  % Could also be obtained by reading from the file name
       if (ibe_pfile<=iminloc) && (iend_pfile>=iminloc)
         n_lower = n;
         ioffset_lower = (iminloc - ibe_pfile)/iskip_DNS;
         inum_lower = min( (iend_pfile - iminloc)/iskip_DNS + 1, iavelen);
       end
       
       if (ibe_pfile<=imaxloc)&&(iend_pfile>=imaxloc)
           n_upper = n;
           if(n_upper>n_lower) 
               inum_upper = (imaxloc - ibe_pfile)/iskip_DNS + 1;
           end
           break
       end
    end
    
    nfile = n_upper - n_lower + 1; % Total number of files that need to be read
    
    
    nn = 0;
    for n=n_lower:n_upper
       nn = nn + 1;
       ibe_pfile  = ibe_DNS + ((n-1)*ibuffer)*iskip_DNS;
       iend_pfile = min(ibe_pfile + (ibuffer-1)*iskip_DNS, iend_DNS);
       fname_irange(nn,:) = sprintf(strcat('i%4.4d-%4.4d'),ibe_pfile,iend_pfile);
       %sprintf(strcat('i','%4.4d-%4.4d'),ibe_pfile,iend_pfile)
    end
    
    TSkplane(1:nfile).gname = '/kplane';
    TSkplane(1:nfile).rank  = 4;
    
    % Files that are fully read
    for nn=1:nfile
       TSkplane(nn).dimsf(1) = ntpoint_total;
       TSkplane(nn).dimsf(2) = (jend_ave - jbe_ave)/jskip_DNS + 1;
       TSkplane(nn).dimsf(3) = ibuffer;
       TSkplane(nn).dimsf(4) = 1;
       TSkplane(nn).offset(1) = 0;
       TSkplane(nn).offset(2) = jbe_ave - 1;
       TSkplane(nn).offset(3) = 0;
       TSkplane(nn).offset(4) = 0;
    end
    
    % Files that are partially read in i-dir
    
    TSkplane(1).dimsf(3)  = inum_lower;
    TSkplane(1).offset(3) = ioffset_lower;
    if nfile>1
       TSkplane(nfile).dimsf(3) = inum_upper; 
    end
   
    
    for nn=1:nfile
       TSkplane(nn).dimsm(1:4)  = TSkplane(nn).dimsf(1:4);
       TSkplane(nn).block(1:4)  = TSkplane(nn).dimsf(1:4);
       TSkplane(nn).count(1:4)  = 1;
       TSkplane(nn).stride(1:4) = 1;    
    end
    
    TSkplane.IsHSInitialized = 'true';

    ntpoint_tmp = 0;
    nnn = imints_k;


    % reading variables
    TSkplane.dname = 'p';
    TSkplane.dnum = 1;
    for n=1:nfile
       TSkplane(n).fname = sprintf(strcat(filepath,'timeseries_kplane%4.4d_',fname_irange(n,:),'.h5'),kloc);  
       buffer(ntpoint_tmp+1:ntpoint_tmp+TSkplane(n).dimsf(1),1:TSkplane(n).dimsf(2),nnn:(nnn+TSkplane(n).dimsf(3)-1),1:TSkplane(n).dimsf(4),1:TSkplane(n).dnum)  ...
          = ReadTSHDF5_4D(TSkplane(n));
       nnn = nnn + TSkplane(n).dimsf(3);
    end

else
    
    ibe_DNS   = DNSIndex_k.ibe;
    iend_DNS  = DNSIndex_k.iend;
    iskip_DNS = DNSIndex_k.iskip;
    jbe_DNS   = DNSIndex_k.jbe;
    jend_DNS  = DNSIndex_k.jend;
    jskip_DNS = DNSIndex_k.jskip;
   
    % begin, end and number of spatial points for Averaging
    if mod(ibe_k-ibe_DNS,iskip_DNS)~=0
       ibe_ave = ibe_k+iskip_DNS - mod(ibe_k-ibe_DNS,iskip_DNS);
    else
       ibe_ave = ibe_k;
    end
    nxpoint_k = (iend_k-ibe_ave)/iskip_DNS + 1;
    iend_ave = ibe_ave + (nxpoint_k-1)*iskip_DNS;
    
    if mod(jbe_k-jbe_DNS,jskip_DNS)~=0
       jbe_ave = jbe_k+jskip_DNS - mod(jbe_k-jbe_DNS,jskip_DNS);
    else
       jbe_ave = jbe_k;
    end
    nypoint_k = (jend_k-jbe_ave)/jskip_DNS + 1;
    jend_ave  = jbe_ave + (nypoint_k-1)*jskip_DNS;
    fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf('DNS Output timeseries spatial index range\n');
    fprintf('ibe = %4g, iend = %4g\n',ibe_DNS,iend_DNS);
    fprintf('jbe = %4g, jend = %4g\n',jbe_DNS,jend_DNS);
    fprintf('Actual Spatial Average range\n');
    fprintf('ibe_ave = %4g, iend_ave = %4g\n',ibe_ave,iend_ave);
    fprintf('jbe_ave = %4g, jend_ave = %4g\n',jbe_ave,jend_ave);
    fprintf('Number of spatial points for Averaging\n');
    fprintf('nxpoint_k = %4g, nypoint_k = %4g\n',nxpoint_k,nypoint_k);
    
    % Dimension range for TSData
    imints_k = 1 - min(ixwindowl, (ibe_ave-ibe_DNS)/iskip_DNS);
    imaxts_k = nxpoint_k + min(ixwindowr, (iend_DNS-iend_ave)/iskip_DNS);
    iavelen = imaxts_k - imints_k + 1;
    fprintf('streamwise Index range for kplane buffer\n');
    fprintf('imints_k = %4g, imaxts_k = %4g, iavelen (buffer length) = %4g\n', imints_k,imaxts_k,iavelen);
    
    iminloc = max(ibe_ave  - ixwindowl*iskip_DNS,ibe_DNS);
    imaxloc = min(iend_ave + ixwindowr*iskip_DNS,iend_DNS);
    fprintf('Physical i-index range read TS (including xwindow)\n');
    fprintf('iminloc = %4g, imaxloc = %4g\n',iminloc,imaxloc);
    fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    
    TSkplane.gname = '/kplane';
    TSkplane.rank  = 4;
    
    TSkplane.dimsf(1) = DNSIndex_k.ntpoint;
    TSkplane.dimsf(2) = (jend_ave - jbe_ave)/jskip_DNS + 1;
    TSkplane.dimsf(3) = imaxts_k - imints_k + 1;
    TSkplane.dimsf(4) = 1;
    TSkplane.offset(1) = 0;
    TSkplane.offset(2) = jbe_ave - jbe_DNS;
    TSkplane.offset(3) = ibe_ave - ibe_DNS - ixwindowl;
    TSkplane.offset(4) = 0;
    
    TSkplane.dimsm(1:4)  = TSkplane.dimsf(1:4);
    TSkplane.block(1:4)  = TSkplane.dimsf(1:4);
    TSkplane.count(1:4)  = 1;
    TSkplane.stride(1:4) = 1; 
    
    TSkplane.IsHSInitialized = 'true';
    
    ntpoint_tmp = 0;
    nnn = imints_k;

    % reading variables
    TSkplane.dname = 'p';
    TSkplane.dnum = 1;
    for n=1:num_file
       TSkplane.fname = sprintf(strcat(filepath,'timeseries_%8.8d.h5'),file_be + (n-1)*file_skip);
       buffer(ntpoint_tmp+1:ntpoint_tmp+TSkplane.dimsf(1),1:TSkplane.dimsf(2),nnn:(nnn+TSkplane.dimsf(3)-1),1:TSkplane.dimsf(4),1:TSkplane.dnum)  ...
          = ReadTSHDF5_4D(TSkplane);
       ntpoint_tmp = ntpoint_tmp + TSkplane.dimsf(1);
    end
    
end % end ireadVOL.eq.0