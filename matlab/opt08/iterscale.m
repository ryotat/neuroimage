function [xx,status]=iterscale(xx, problem, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'epsf', 1e-6, ...
                      'epsg', 1e-5, ...
                      'epssvd', 1e-14, ...
                      'maxiter', 0,...
                      'max_linesearch', 30,...
                      'nmexp',0,...
                      'logfile', [],...
                      'display', 1,...
                      'D0', []);

if ~isempty(opt.logfile)
  problem.fid = fopen(opt.logfile, 'w');
end

problem.nfeval = 0;
problem.t0     = cputime;

D = feval(problem.init_scaling, problem, opt);

%% Generate differentiable subproblem
prob = packvars(problem, problem.Y, problem.lambda, 'squared', problem.loss(5:end), 'fro');

kk = 1;

while 1
  %% Solve a differentiable subproblem
  prob.X      = feval(problem.apply_scaling, D, problem);
  [xx,stat]=feval(prob.solver, xx, prob, opt, 'epsf',1e-10,'maxiter',0,'display', 0);
  
  problem.nfeval = problem.nfeval + stat.nfeval;

  xx = feval(problem.apply_scaling, D, xx, problem);

  %% Compute the scaling coefficient
  D = feval(problem.adjust_scaling, xx, problem);

  if ~isreal(D)
    keyboard;
  end
  
  %% Compute the objective
  problem.lambda = prob.lambda/feval(problem.norm, xx, problem)^opt.nmexp;
  [fval,gg,dval,problem]=feval(problem.obj, xx, problem, opt);

  
  % Progress report
  feval(problem.progress, kk, xx, fval, gg, dval, problem.dpen, problem, 0, 0, 0, opt);

  if fval-dval<opt.epsf*max(max(fval,dval),1) || max(abs(gg))<opt.epsg
    fprintf('Optimization success! gap=%g\n', fval-dval);
    ret=0;
    break;
  end
  
  if kk==opt.maxiter
    fprintf('Maximum #iterations=%d reached.\n', kk);
    ret = -3;
    break;
  end

  kk = kk + 1;
end


if ~isempty(opt.logfile)
  fclose(problem.fid);
end


status=struct('ret', ret,...
              'kk', kk,...
              'fval', fval,...
              'dval', dval,...
              'gg', gg,...
              'nfeval', problem.nfeval,...
              'time', cputime-problem.t0,...
              'opt', opt);
