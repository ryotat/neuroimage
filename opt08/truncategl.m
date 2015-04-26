function ww = truncategl(ww, problem, opt)
  
if problem.hasbias
  bias = ww(end);
  ww=ww(1:end-1);
else
  bias = [];
end

ww=reshape(ww,[problem.ns, problem.nc]);

nm  = sqrt(sum(ww.^2));
ix0 = find(nm<=opt.epsh);

if length(ix0)>0
  fprintf('Zeroing %d components:', length(ix0));
end


ww(:,ix0)=0;


ww = [reshape(ww,[problem.ns*problem.nc,1]); bias];
