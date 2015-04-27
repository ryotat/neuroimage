function dval = dnormds(ww, gg, problem, opt)
if problem.hasbias
  gg = gg(1:end-1);
end

lennm = min([problem.ns; problem.nc]);
ss = zeros([1, sum(lennm)]);
C  = zeros([1, sum(lennm)]);

ix0=0;
ixnm0=0;
for ii=1:length(lennm)
  sz  = [problem.ns(ii), problem.nc(ii)];
  ix  = ix0+(1:prod(sz));
  ix0 = ix(end);

  ixnm=ixnm0+(1:lennm(ii));
  ixnm0=ixnm(end);
  ss(ixnm)=svd(reshape(gg(ix),sz))';
  C(ixnm) =problem.C(ii)*ones(1,lennm(ii));
end

dval = -max(ss./C);
  