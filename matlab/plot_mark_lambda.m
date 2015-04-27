function ix=plot_mark_lambda(xvalue, lmd, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'logscale',1,'lambda', xvalue);

if opt.logscale
  xx = log(xvalue);
else
  xx = xvalue;
end


ix = findnearest(opt.lambda, lmd);
if ~isempty(ix)
  plot(ones(2,1)*xx(ix), ylim'*ones(1,length(ix)),'--', 'color',[.5 .5 .5],'linewidth',2);
end

fprintf('lmd=%g ix=%d nm=%g\n', lmd, ix, xvalue(ix));
