
% initial timeseries volume info

function fsol = InitFlowHDF52D()
  fsol.gname = '/';
  fsol.sname = 'time';
  fsol.rank = 2;
  fsol.count = [1 1];
  fsol.stride = [1 1];
  fsol.offset = [0 0];
end
