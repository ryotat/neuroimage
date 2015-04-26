function D = adjust_scaling_l1_wip(xx, prob)

dd = prob.ns*prob.nc;

if prob.hasbias
  xx = xx(1:end-1);
end

D = sqrt(abs(xx));

