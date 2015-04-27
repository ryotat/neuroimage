function [fv, Ws, Wt] = proc_whitenST(fv, method, varargin)
% proc_whitenST - whitening both space and time
%
% [epo, Ws, Wt] = proc_whitenST(epo)
%

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'exps', -0.5, 'expt', -0.5);

if ~exist('method','var') || isempty(method)
  method = 'st';
end

if iscell(method)
  opt.exps = method{2};
  opt.expt = method{3};
  method = method{1};
end

[T,C,n]=size(fv.x);

if strcmp(method, 'st')
  fprintf('whitening: space&time [exps=%g, expt=%g]\n',opt.exps,opt.expt);
  SigmaS = covPooled(fv.x);
  SigmaT = covPooled(permute(fv.x, [2,1,3]));

  if opt.exps==-0.5 && opt.expt==-0.5
    Ws=inv(sqrtm(SigmaS));
    Wt=inv(sqrtm(SigmaT));
  else
    Ws=SigmaS^opt.exps;
    Wt=SigmaT^opt.expt;
  end
  
  fv.x=prodAt(fv.x, Wt, 1);
  fv.x=prodAt(fv.x, Ws, 2);
elseif strcmp(method, 's')
  fprintf('whitening: space only\n');
  [fv, Ws] = proc_whiten(fv);
  Wt = eye(T);
elseif strcmp(method, 't')
  fprintf('whitening: time only\n');
  SigmaT = covPooled(permute(fv.x, [2,1,3]));
  Ws = eye(C);
  Wt = inv(sqrtm(SigmaT));
  fv.x = prodAt(fv.x, Wt, 1);
elseif strcmp(method, 'scalest')
  fprintf('whitening: space&time scaling only\n');
  vs = sqrt(mean(std(fv.x).^2,3));
  vt = sqrt(mean(std(fv.x,[],2).^2,3));
  Ws = diag(1./vs);
  Wt = diag(1./vt);
  fv.x = prodAt(fv.x, Wt, 1);
  fv.x = prodAt(fv.x, Ws, 2);
elseif strcmp(method, 'scales+t')
  fprintf('whitening: scale space & whiten time\n');
  vs = sqrt(mean(std(fv.x).^2,3));
  SigmaT = covPooled(permute(fv.x, [2,1,3]));
  Ws = diag(1./vs);
  Wt = inv(sqrtm(SigmaT));
  fv.x = prodAt(fv.x, Wt, 1);
  fv.x = prodAt(fv.x, Ws, 2);
elseif strcmp(method, 'scales+fourier') && strcmp(opt.group,'fourier')
  fprintf('whitening: scale space & time fourier\n');
  vs = sqrt(mean(std(fv.x).^2,3));
  Ws = diag(1./vs);
  Wt = fft(eye(T))';
  fv.x = prodAt(fv.x, Wt', 1);
  fv.x = prodAt(fv.x, Ws, 2);
elseif strcmp(method, 'concats')
  fprintf('whitening: spcae (concat)\n');
  [fv, Ws] = proc_whiten_concat(fv);
  Wt = eye(T);
else
  fprintf('whitening: none\n');
  Ws = eye(C);
  Wt = eye(T);
end
