function ww = truncateds(ww, problem, opt)
  
if problem.hasbias
  bias = ww(end);
  ww=ww(1:end-1);
else
  bias = [];
end

[U,S,V]=svd(reshape(ww,[problem.ns, problem.nc]));

nm  = diag(S);
M   = max(nm);
nmp = nm.*(nm>M*opt.epsh);

fprintf('Zeroing %d components:', sum(nm<=M*opt.epsh));


if problem.ns<problem.nc
  ww = U*spdiag(nmp)*V(:,1:problem.ns)';
else
  ww = U(:,1:problem.nc)*spdiag(nmp)*V';
end

ww = [reshape(ww,[problem.ns*problem.nc,1]); bias];
