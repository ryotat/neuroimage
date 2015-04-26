function D = spdiag(d)

if isempty(d)
  D = [];
  return;
end

if size(d,1)<size(d,2)
  d = d';
end

D = spdiags(d,0,size(d,1),size(d,1));