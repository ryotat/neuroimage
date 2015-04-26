function D = init_scaling_l1_wip(prob, opt)

dd = prob.ns*prob.nc;

if ~isempty(opt.D0);
  D = opt.D0;
else
  D = ones(dd,1);
end
