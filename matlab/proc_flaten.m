function dat= proc_flaten(dat)
%dat= proc_flaten(dat)
%
% reshape data matrix to data vector (clash all but last dimensions)

% bb, ida.first.fhg.de


if isnumeric(dat),
  % Also flaten the data if it is a plain data matrix
  sz = size(dat);
  dat = reshape(dat, [prod(sz(1:end-1)) sz(end)]); 
elseif isstruct(dat),
  % Old code from the BCI toolbox:
  if isstruct(dat.x),
    dat= proc_flatenGuido(dat);
  else
    sz = size(dat.x);
    dat.x = reshape(dat.x, [prod(sz(1:end-1)) sz(end)]);
  end
else
  error('Don''t know how to flatten this type of data');
end



% Dec 2008: copied from IDA toolbox. 
% All rights belong to the authors and Fraunhofer FIRST.IDA.
% http://ida.first.fraunhofer.de/homepages/ida/


