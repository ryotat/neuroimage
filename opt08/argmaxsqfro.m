function dder = argmaxsqfro(dd, ww, fval, gg, problem, opt)

dder = dd'*gg;

if problem.hasbias
  ww = ww(1:end-1);
  dd = dd(1:end-1);
end

ix0=0;
for ii=1:length(problem.lambda)
  sz  = [problem.ns(ii), problem.nc(ii)];
  ix  = ix0+(1:prod(sz));
  ix0 = ix(end);
  
  dder = dder + problem.lambda(ii)*dd(ix)'*ww(ix);
end


