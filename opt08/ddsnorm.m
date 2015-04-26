function dval = ddsnorm(ww, gg, problem, opt)
if problem.hasbias
  gg = gg(1:end-1);
end

ss = svd(reshape(gg, [problem.ns, problem.nc]));
dval = -max(ss'./problem.C);
  