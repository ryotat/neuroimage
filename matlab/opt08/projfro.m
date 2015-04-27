function ww = projfro(ww, problem, opt)

if problem.hasbias
  bias = ww(end);
  ww=ww(1:end-1);
else
  bias = [];
end

len = problem.ns.*problem.nc;

ww = shrink_l2weight(ww, problem.C, len, opt.display-100);

ww = [ww; bias];
