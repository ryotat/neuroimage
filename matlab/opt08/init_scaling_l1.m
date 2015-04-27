function D = init_scaling_l1(prob, opt)

dd = prob.ns*prob.nc;

if ~isempty(opt.D0);
  D = opt.D0;
else
  D = sqrt(ones(dd,1)/dd);
end
