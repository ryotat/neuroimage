function dval = dglnorm(ww, gg, problem, opt)
  
dval = -max(sqrt(sum(reshape(gg,[problem.ns,problem.nc]).^2))./problem.C);
  