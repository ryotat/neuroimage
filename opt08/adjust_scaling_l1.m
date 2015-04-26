function D = adjust_scaling_l1(xx, prob)

dd = prob.ns*prob.nc;

if prob.hasbias
  xx = xx(1:end-1);
end

xabs = abs(xx);
D = sqrt(xabs/sum(xabs));

