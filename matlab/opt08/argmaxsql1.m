function [dder,problem] = argmaxsql1(dd, ww, fval, gg, problem, opt)

dder = dd'*gg;

if problem.hasbias
  ww = ww(1:end-1);
  dd = dd(1:end-1);
end

nm = abs(ww);

lmnm = problem.lambda*sum(nm);

ix1 = find(nm>0);
ix0 = find(nm==0);

dder = dder + lmnm*sum(dd(ix1).*(ww(ix1)./nm(ix1)));

if ~isempty(ix0)
  dder = dder + lmnm*sum(abs(dd(ix0)));
end

