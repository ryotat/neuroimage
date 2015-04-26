% train_gen - Trains general models for BCI with different
%             optimization routines
%
% Syntax:
%  C = train_gen(xTr, yTr, lambda, <opt>)
%
% Input:
%  xTr   - Input EEG [nSamples,nChannels,n] or [nSamples*nChannels,n]
%  yTr   - Labels    [ncls,n] or [2,ncls*n]
% lambda - Regularization constant
%  opt   - Struct or property/value list of optional properties:
%   .nSamples  - number of sampled time-points
%   .nChannels - number of channels
%   .ncls      - number of classes    (default 2)
%   .whitening - whitening option see proc_whitenST (default 'st')
%   .loss      - loss function        (default 'lr')
%   .reg       - regularizer          (default 'ds')
%   .group     - group dimension (space/time/fourier)
%                                     (default '')
%   .solver    - optimization routine (default 'projgrad')
%   .display   - display level        (default 'final')
%   .W0        - initial value        (default [])
%   .epsf      - relative duality gap tolerance (default 1e-4)
%   .maxiter   - maximum number of iterations   (default 100000)
%   (any other option similarly applies to the optimization routines)
% Output:
%  C     - Classifier
%
% See also:
%  trainClassifier proc_whitenST packvars projgrad sublbfgs
%
% Ryota Tomioka 2009

function cls = train_gen(xTr, yTr, lambda, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'whitening', 'st',...
                      'loss','lr',...
                      'reg', 'ds',...
                      'group','',...
                      'solver', 'projgrad',...
                      'ncls', 2,...
                      'W0', [],...
                      'display', 'final',...
                      'epsf',1e-4,...
                      'maxiter',100000);

if isempty(lambda)
  error('Lambda must be specified.');
end

if ~isnumeric(opt.display)
  opt.display = find(strcmp({'none','final','iter'}, opt.display))-1;
end


xTr = shiftdim(xTr);

sz = size(xTr);
if sz(end)~=size(yTr,2)
  error('Number of samples mismatch!');
end

iValid = find(any(yTr));

if ndims(xTr)<3
  T = opt.nSamples;
  C = opt.nChannels;
  n = sz(end);
  xTr = reshape(xTr, [T, C, n]);
else
  [T,C,n]=size(xTr);
end

fv = struct('x', xTr(:,:,iValid), 'y', yTr(:,iValid));

if size(yTr,1)~=opt.ncls
  nTrials = length(iValid)/opt.ncls;
else
  nTrials = length(iValid);
end

%% Whitening
[fv,Ws,Wt] = proc_whitenST(fv, opt.whitening);

if ~isequal(opt.reg,'gl') && ~isempty(opt.group)
  error('Regularization must be ''gl'' to specify grouping.')
end

if strcmp(opt.group,'time')
  fprintf('grouping: time\n');
  xTr = permute(fv.x, [2,1,3]);
elseif strcmp(opt.group,'fourier')
  fprintf('grouping: fourier\n');
  xTr = reshape(fv.x, [1, T*C, length(iValid)]);
  xTr = [real(xTr); imag(xTr)];
else
  if isequal(opt.reg,'gl')
    fprintf('grouping: space\n');
  end
  xTr = fv.x;
end

switch(opt.solver)
case  'projgrad'
  algo = 'constraint';
case 'sublbfgs'
  algo = 'squared';
 case 'iterscale'
  algo = 'iterscale';
 otherwise
  algo = opt.solver;
end

%% Initial value
if ~isempty(opt.W0)
  if size(opt.W0)==[T,C]
    xx0=reshape(opt.W0,[T*C,1]);
  else
    xx0=opt.W0;
  end
else
  xx0=zeros(T*C,1);
end

%% Call the solver

fprintf('loss=%s reg=%s\n',opt.loss, opt.reg);
if size(fv.y,1)~=opt.ncls
  Y=reshape(fv.y(2,:), [opt.ncls, nTrials]);
else
  Y=2*fv.y(2,:)-1;
end

problem = packvars(fv.x, Y, lambda, algo, opt.loss, opt.reg)

if problem.hasbias && length(xx0)<size(problem.X,1)
  xx0 = [xx0;0];
end

[xx, status] = feval(problem.solver, xx0, problem, opt);

if problem.hasbias
  W = reshape(xx(1:end-1),[T,C]);
else
  W = reshape(xx,[T,C]);
end

if strcmp(opt.group,'time')
  W=W';
elseif strcmp(opt.group,'fourier')
  Wr = reshape(W(1,:), [T,C]); Wi = reshape(W(2,:), [T,C]);
  W=Wr+sqrt(-1)*Wi;
end

if problem.hasbias
  cls = struct('W',W,'bias',xx(end),'Ws',Ws,'Wt',Wt,'status',status);
else
  cls = struct('W',W,'Ws',Ws,'Wt',Wt,'status',status);
end



fprintf('Elapsed time=%g\n', status.time);





