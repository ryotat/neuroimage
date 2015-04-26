function m = medin(x)

x=sort(x);

if size(x,1)==1
  ix=ceil(length(x)/2);
  m = x(ix);
else
  ix=ceil(size(x,1)/2);
  m = x(ix,:);
end

