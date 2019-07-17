

function DNSIndex_k = ReadDNS_index_kplane(filepath)
  % Input: filepath
  
  % Output: Structure array
  fname = strcat(filepath,'series_time_ascii.dat');
  fID=fopen(fname,'rt'); 
  for n=1:2
    time(n) = fscanf(fID,'%e',1);
  end
  fclose(fID);
  DNSIndex_k.dt = time(2) - time(1);  % dt_sample
  
  DNSIndex.gname = '/kplane';
  Dname = {'dx' 'dy' 'istart' 'iend' 'iskip' 'jstart' 'jend' 'jskip' 'nkplane' 'ntpoint'};
  DNSIndex.fname = strcat(filepath,'DNS_index.h5');
  for n=1:10
     DNSIndex.sname = Dname{n};
     dset(n) = ReadScalar(DNSIndex);
  end
  DNSIndex_k.dx = dset(1);
  DNSIndex_k.dy = dset(2);
  DNSIndex_k.ibe = dset(3);
  DNSIndex_k.iend = dset(4);
  DNSIndex_k.iskip = dset(5);
  DNSIndex_k.jbe = dset(6);
  DNSIndex_k.jend = dset(7);
  DNSIndex_k.jskip = dset(8);
  DNSIndex_k.nkplane = dset(9);
  DNSIndex_k.ntpoint = dset(10);
  
  DNSIndex.rank = 1;
  DNSIndex.count = ones(1,DNSIndex.rank);
  DNSIndex.stride = ones(1,DNSIndex.rank);
  DNSIndex.offset = zeros(1,DNSIndex.rank);
  DNSIndex.dimsm = [DNSIndex_k.nkplane];
  
  DNSIndex.dname = 'klocs';
  DNSIndex_k.klocs = ReadHDF5(DNSIndex);
   
end 