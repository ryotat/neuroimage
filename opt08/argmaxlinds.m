function dder = argmaxlinds(dd, ww, fval, gg, problem, opt)

dder = dd'*gg;

if problem.hasbias
  ww = ww(1:end-1);
  P = reshape(dd(1:end-1), [problem.ns, problem.nc]);
else
  P = reshape(dd, [problem.ns, problem.nc]);
end
ww = reshape(ww, [problem.ns, problem.nc]);

[U,S,V] = svd(ww);
nm      = diag(S);

ix1 = find(nm>opt.epssvd);
ix0 = find(nm<=opt.epssvd);

U1=U(:,ix1); V1=V(:,ix1);
U0=U(:,ix0); V0=V(:,ix0);


dder = dder + problem.lambda*trace(U1'*P*V1);

if ~isempty(ix0)
  dder = dder + problem.lambda*sum(svd(U0'*P*V0));
end

