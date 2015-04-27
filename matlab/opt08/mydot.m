function d = mydot(x, y)
  
d=0;  
for ii=1:length(x)
  d = d + x(ii)*y(ii);
end
