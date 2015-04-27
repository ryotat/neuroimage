function dder = argmaxsqgl(dd, ww, fval, gg, problem, opt)

dder = dd'*gg;

if problem.hasbias
  ww = ww(1:end-1);
  P = reshape(dd(1:end-1), [problem.ns, problem.nc]);
else
  P = reshape(dd, [problem.ns, problem.nc]);
end
ww = reshape(ww, [problem.ns, problem.nc]);

nm = sqrt(sum(ww.^2));

lmnm = problem.lambda*sum(nm);

ix1 = find(nm>0);
ix0 = find(nm==0);

dder = dder + lmnm*sum(sum(P(:,ix1).*(ww(:,ix1)/spdiag(nm(ix1)))));

if ~isempty(ix0)
  dder = dder + lmnm*sum(sqrt(sum(P(:,ix0).^2)));
  % dder2 = dder + lmnm*sum(norm(P(:,ix0)));;
  % fprintf('\ndder1=%g dder2=%g\n',dder1,dder2);
end

