% sublbfgs - Subgradient based L-BFGS algorithm
%
% Syntax:
%  [xx, status] = sublbfgs(xx, problem, <opt>)
%
% Input:
%  xx      - Initial point for optimization
%  problem - Optimization problem with fields:
%   .obj        - objective function
%   .projection - truncate singular values smaller than epsh
%   .argmax     - calculate directional derivative
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

function [xx, status] = sublbfgs(xx, problem, varargin)
  
opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'm', 6,...
                      'epsh', 1e-6, ...
                      'epsf', 1e-6, ...
                      'epsg', 1e-5, ...
                      'epssvd', 1e-14, ...
                      'ftol', 1e-5, ...
                      'maxiter', 0,...
                      'max_linesearch', 30,...
                      'logfile', [],...
                      'display', 1,...
                      'maxdval',0);

if ~isempty(opt.logfile)
  problem.fid = fopen(opt.logfile, 'w');
end

nn = size(xx,1);

if ~isfield(problem,'dpen')
  problem.dpen=[];
end

problem.nfeval = 0;
problem.t0 = cputime;

% Limited memory
lm = repmat(struct('s',zeros(nn,1),'y',zeros(nn,1),'ys',0,'alpha',0),[1, opt.m]);

[fval,gg,dval,problem,floss,gloss]=feval(problem.obj, xx, problem, opt);

maxdval = dval;

% The initial step is gradient
dd = -gg;

kk = 1;
stp = 1/norm(dd);
step_grad = 1.0;

bResetLBFGS = 0;
ixend = 1;
bound = 0;
while 1
  fp = fval;
  xxp = xx;
  ggp = gg;

  % Projection
  dd=feval(problem.projection, xx+step_grad*dd, problem, opt)-xx;

  
  dder = feval(problem.argmax, dd, xx, floss, gloss, problem, opt);
  % fprintf('dder=%.20f:', dder);
  if dder>0
    dd = -gg;
    bResetLBFGS = 1;
    dder = -norm(gg);
  end

  % Perform line search
  [ret, xx,fval,gg,dval,problem,stp,floss,gloss]=...
      feval(problem.linesearch, xx, fval, gg, dd, dder, stp, problem, opt);
  
  if opt.maxdval
    maxdval = max(maxdval, dval);
  else
    maxdval = dval;
  end
  
  
  if ret<0
    break;
  end

  
  
  % Evaluate again
  % [fval,gg,dval,problem]=feval(problem.obj, xx, problem, opt);

  
  % Progress report
  feval(problem.progress, kk, xx, fval, gg, dval, problem.dpen, problem, bound, stp, step_grad, ...
        opt);

  
% $$$   if stp<1.0
% $$$     step_grad = step_grad/2;
% $$$     fprintf('step_grad=%g\n', step_grad);
% $$$   else
% $$$     step_grad = step_grad*sqrt(2);
% $$$   end

  if fval-maxdval<opt.epsf*max(max(fval,maxdval),1) || max(abs(gg))<opt.epsg
    if opt.display>0
      fprintf('Optimization success! gap=%g\n', fval-maxdval);
    end
    ret=0;
    break;
  end
  
%  if fp-fval<opt.epsf*max(max(fval,fp),1)
%     bResetLBFGS = 1;
%   end
  
  if kk==opt.maxiter
    if opt.display>0
      fprintf('Maximum #iterations=%d reached.\n', kk);
    end
    ret = -3;
    break;
  end


  % if abs(fval-fp)<1e-9*max(max(fval,fp),1)
  % fprintf('Resetting L-BFGS matrix.\n');
  %    bResetLBFGS = 1;
  % end
  
  if bResetLBFGS
    fprintf('Resetting LBFGS matrix.\n');
    bound = 0;
    ixend = 1;
    bResetLBFGS = 0;
    ys = 1; yy = 1;
  else
    % L-BFGS update
    if opt.m>0
      lm(ixend).s = xx-xxp;
      lm(ixend).y = gg-ggp;
      ys = lm(ixend).y'*lm(ixend).s; yy = sum(lm(ixend).y.^2);
      lm(ixend).ys  = ys;
    else
      ys = 1; yy = 1;
    end
    
    bound = min(bound+1, opt.m);
    ixend = (opt.m>0)*(mod(ixend, opt.m)+1);
  end
  
  % Initially set the negative gradient as descent direction
  dd = -gg;
  
  jj = ixend;
  for ii=1:bound
    jj = mod(jj + opt.m -2, opt.m)+1;
    lm(jj).alpha = lm(jj).s'*dd/lm(jj).ys;
    dd = dd -lm(jj).alpha*lm(jj).y;
  end

  dd = dd *(ys/yy);
  
  for ii=1:bound
    beta = lm(jj).y'*dd/lm(jj).ys;
    dd = dd + (lm(jj).alpha-beta)*lm(jj).s;
    jj = mod(jj,opt.m)+1;
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
              'opt', opt);
