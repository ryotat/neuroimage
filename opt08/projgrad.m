% projgrad - Gradient projection method for constrained
%            optimization
% Syntax:
%  [xx, status] = projgrad(xx, problem, varargin)
%
% Input:
%  xx      - Initial point for optimization
%  problem - Optimization problem with fields:
%   .obj        - objective function
%   .dpen       - dual penalty function
%   .projection - projection to the feasible set
%   .linesearch - line search routine
%   .progress   - report progress
%  opt     - Struct or property/value list of optional properties:
%   .m          - size of limited memory
%   .epsh       - threshold for small singular values
%   .epsf       - relative duality gap tolerance
%   .epsg       - gradient tolerance
%   .maxiter    - maximum number of iterations
%   .logfile    - name of logfile
%   .display    - display level
% 
% Output:
%  xx      - Final point of optimization
%  status  - Status (e.g., primal/dual objective value,
%            gradient,...)
%
% See also:
%  packvars
%
% Ryota Tomioka 2009
function [xx, status] = projgrad(xx, problem, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'epsf', 1e-6, ...
                      'ftol', 1e-12, ...
                      'maxiter', 0,...
                      'max_linesearch', 20,...
                      'logfile', [],...
                      'display', 1,...
                      'maxdval',0);

if ~isempty(opt.logfile)
  problem.fid = fopen(opt.logfile, 'w');
end

problem.nfeval = 0;
problem.t0 = cputime;

[fval,gg,dval,problem]=feval(problem.obj, xx, problem, opt);
dpen = feval(problem.dpen, xx, gg, problem, opt);
dval = dval + dpen;

maxdval = dval;

kk = 1;
stp = 1/norm(gg);
step_grad = stp;
while 1
  fp = fval;
  
  % The step is gradient
  dd = -gg;

  % Take a step and project
  sg = step_grad;
  cc=0;
  while 1
    dd = feval(problem.projection, xx+sg*dd, problem, opt) - xx;

    if gg'*dd<0
      break;
    end
    if cc>opt.max_linesearch
      fprintf('cc=%d sg=%g Cannot find descent direction\n',cc,sg);
      break;
    end
    sg = sg/2;
    cc = cc + 1;
  end

  % Perform line search
  stp = 1.0;
  [ret, xx,fval,gg,dval,problem,stp]=...
      feval(problem.linesearch, xx, fval, gg, dd, gg'*dd, stp, problem, opt);
  dpen = feval(problem.dpen, xx, gg, problem, opt);
  dval = dval + dpen;

  if opt.maxdval
    maxdval = max(maxdval,dval);
  else
    maxdval = dval;
  end

  
  if stp<1.0
    step_grad = step_grad/2;
    % fprintf('step_grad=%g\n', step_grad);
  else
    step_grad = step_grad*2^(min(0.5,abs(fval-dval)));
  end
  
  % Progress report
  feval(problem.progress, kk, xx, fval, gg, dval, dpen, problem, 0, stp, step_grad, opt);

  if fval<dval
    keyboard;
end

  
  if fval>=maxdval && fval-maxdval<opt.epsf*max(max(fval,maxdval),1);
    fprintf('Optimization success! gap=%g\n', fval-maxdval);
    ret=0;
    break;
  end
  
  if kk==opt.maxiter
    fprintf('Maximum #iterations=%d reached. fval=%g dval=%g\n', kk, ...
            fval, maxdval);
    ret = -3;
    break;
  end

  if ret<0
    break;
  end
  

  stp = 1.0;
  
  kk = kk + 1;
end


if ~isempty(opt.logfile)
  fclose(problem.fid);
end


status=struct('ret', ret,...
              'kk', kk,...
              'fval', fval,...
              'dval', maxdval,...
              'gg', gg,...
              'nfeval', problem.nfeval,...
              'time', cputime-problem.t0,...
              'opt',opt);
