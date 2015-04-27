function ww = projgl(ww, problem, opt)
ww = reshape(ww,[problem.ns, problem.nc]);  
nm = sqrt(sum(ww.^2));

nmp = projection_base(nm, problem.C);

ix = nmp>0;

ww(:,ix)  = ww(:,ix)*spdiag((nmp(ix)./nm(ix))');
ww(:,~ix) = 0;

ww = reshape(ww,[problem.ns*problem.nc,1]);