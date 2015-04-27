function [b,a]=getbutter(order, band, fs)

try
  if length(band)==2
    file = sprintf('butter_%d_%d_%d.mat', band(1), band(2), fs);
  else
    file = sprintf('butter_%d_%d.mat', band, fs);      
  end
  load(file);
catch
  [b,a]=butter(order, band/fs*2);
end

