function ww = projds(ww, problem, opt)

if problem.hasbias
  bias = ww(end);
  ww=ww(1:end-1);
else
  bias = [];
end

lennm = min([problem.ns; problem.nc]);
nm = zeros([1, sum(lennm)]);
C  = zeros([1, sum(lennm)]);
U = cell(1,length(lennm));
V = cell(1,length(lennm));

ix0=0;
ixnm0=0;
for ii=1:length(lennm)
  sz  = [problem.ns(ii), problem.nc(ii)];
  ix  = ix0+(1:prod(sz));
  ix0 = ix(end);

  [U{ii},S,V{ii}]=svd(reshape(ww(ix),sz));
  ixnm=ixnm0+(1:lennm(ii));
  ixnm0=ixnm(end);
  nm(ixnm)=diag(S)';
  C(ixnm) =problem.C(ii)*ones(1,lennm(ii));
end

nmp = projection_base(nm, C);


ix0=0;
ixnm0=0;
for ii=1:length(lennm)
  sz  = [problem.ns(ii), problem.nc(ii)];
  ix  = ix0+(1:prod(sz));
  ix0 = ix(end);
  ixnm=ixnm0+(1:lennm(ii));
  ixnm0=ixnm(end);

  if sz(1)<sz(2)
    ww(ix) = U{ii}*spdiag(nmp(ixnm))*V{ii}(:,1:sz(1))';
  else
    ww(ix) = U{ii}(:,1:sz(2))*spdiag(nmp(ixnm))*V{ii}';
  end
end


ww = [ww; bias];
