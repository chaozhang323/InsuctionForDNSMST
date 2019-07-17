
% initial timeseries volume info

function fsol = InitFlowHDF5(ndim)
  fsol.gname = '/';
  fsol.sname = 'time';
  fsol.rank = ndim;
  fsol.count = ones(1,ndim);
  fsol.stride = ones(1,ndim);
  fsol.offset = zeros(1,ndim);
end
