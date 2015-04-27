function dval = dnorml1(ww, gg, problem, opt)

dval = -max(abs(gg)./problem.C);
  