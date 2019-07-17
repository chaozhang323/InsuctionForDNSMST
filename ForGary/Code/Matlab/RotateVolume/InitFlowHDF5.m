
% initial timeseries volume info

function fsol = InitFlowHDF5()
  fsol.gname = '/';
  fsol.sname = 'time';
  fsol.rank = 3;
  fsol.count = [1 1 1];
  fsol.stride = [1 1 1];
  fsol.offset = [0 0 0];
end
