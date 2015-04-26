% packvars - Defines an optimization problem
%
% Syntax:
%  problem = packvars(X, Y, lambda, algo, loss, reg)
%
% Input:
%  X       - Input EEG [nSamples, nChannels, n]
%  Y       - Labels [ncls,n] with zero/one elements or [1,n]
%            with elements from {0,1,...,ncls-1}
%  lambda  - Regularization constant
%  algo    - Algorithm:
%    squared    - squared penalty formulation
%    constraint - norm constraint formulation
%    linear     - linear penalty formulation
%    iterscale  - scaled Frobenius formulation
%  loss    - Loss function:
%    lr         - binary logistic loss
%    mn         - multinomial loss for P300
%    sq         - squared loss
%  reg     - Regularization:
%    l1         - lasso
%    gl         - group lasso
%    ds         - dual spectral (trace norm)
% Output:
%  problem - optimization problem
%
% See also:
%  train_gen sublbfgs projgrad
%
% Ryota Tomioka 2009

function problem = packvars(X, Y, lambda, algo, loss, reg, regopt)

if ~exist('regopt','var')
  regopt = '';
end

if iscell(X)
  [problem.X,ns,nc,nall]=pack_matrices(X{:});
  
elseif isstruct(X)
  problem.X = X.X;
  ns        = X.ns;
  nc        = X.nc;
  nall      = size(X.X,2);
else
  [ns, nc, nall]= size(X);
  problem.X = reshape(X,[ns*nc,nall]);
end

  
n             = size(Y,2);
  
problem.ns    = ns;
problem.nc    = nc;
problem.ncls  = nall/n;



switch(algo)
 case 'constraint'
  problem.C = lambda;
  problem.obj        = ['loss' loss];
  problem.dpen       = ['dnorm' reg];
  problem.projection = ['proj' reg];
  problem.progress   = ['progress' reg];
  problem.suf        = [loss '_' reg '_con'];
  problem.solver     = 'projgrad';
 case 'squared'
  problem.lambda     = lambda;
  problem.obj        = ['objsq' reg];
  problem.loss       = ['loss' loss];
  problem.dpen       = [];
  problem.projection = ['truncate' reg];
  problem.progress   = ['progress' reg];
  problem.argmax     = ['argmaxsq' reg];
  problem.suf        = [loss '_' reg '_sub'];
  problem.solver     = 'sublbfgs';
 case 'linear'
  problem.lambda     = lambda;
  problem.obj        = ['objlin' reg];
  problem.loss       = ['loss' loss];
  problem.dpen       = [];
  problem.projection = ['truncate' reg];
  problem.progress   = ['progress' reg];
  problem.argmax     = ['argmaxlin' reg];
  problem.norm       = ['norm' reg];
  problem.suf        = [loss '_' reg '_sub'];
  problem.solver     = 'sublbfgs';
 case 'iterscale'
  problem.lambda         = lambda;
  if isequal(regopt,'wip')
    problem.obj          = ['objlin' reg];
  else
    problem.obj          = ['objsq' reg];
  end
  problem.loss           = ['loss' loss];
  problem.dpen           = [];
  problem.progress       = ['progress' reg];
  problem.apply_scaling  = ['apply_scaling_' reg];
  if ~isempty(regopt)
    problem.init_scaling   = ['init_scaling_' reg '_' regopt];
    problem.adjust_scaling = ['adjust_scaling_' reg '_' regopt];
  else
    problem.init_scaling   = ['init_scaling_' reg];
    problem.adjust_scaling = ['adjust_scaling_' reg];
  end
  
  problem.suf            = [loss '_' reg '_its'];
  problem.solver         = 'iterscale';
end

switch(loss)
 case 'mn'
  if isequal(size(Y),[problem.ncls, n])
    problem.Y = sparse(Y);
  else
    problem.Y = sparse(Y+1, 1:n, ones(1,n), problem.ncls, n);
  end
  
  problem.hasbias = 0;
 case 'lr'
  if isequal(unique(Y),[0,1])
    problem.Y = 2*Y-1;
  else
    problem.Y = Y;
  end
  problem.hasbias = 1;
  problem.ncls = [];
 case 'sq'
  problem.Y = Y;
  problem.hasbias = 0;
end


if problem.hasbias && ~isstruct(X)
  problem.X = [problem.X;ones(1,nall)];
end


problem.linesearch = 'linesearch_backtracking';

problem.klog = [];

