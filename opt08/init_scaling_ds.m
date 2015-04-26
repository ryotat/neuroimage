function D = init_scaling_ds(prob, opt)

if ~isempty(opt.D0);
  D = opt.D0;
else
  D = sqrt(eye(prob.ns)/prob.ns);
end
