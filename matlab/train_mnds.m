function cls = train_mnds(xTr, yTr, lambda, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'whitening', 'st',...
                      'reg', 'ds',...
                      'solver', 'projgrad',...
                      'ncls', 6,...
                      'W0', [],...
                      'display', 'final');

if isempty(lambda)
  error('Lambda must be specified.');
end

% if ~isnumeric(opt.display)
%  opt.display = find(strcmp({'none','final','iter'}, opt.display))-1;
% end


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
nTrials = length(iValid)/opt.ncls;

%% Whitening
[fv,Ws,Wt] = proc_whitenST(fv, opt.whitening);


%% Initial value
if size(opt.W0)==[T,C]
  xx0=reshape(opt.W0,[T*C,1]);
else
  xx0=zeros(T*C,1);
end

switch(opt.solver)
 case 'dal'
  xx0=reshape(xx0,[T,C]);
  X=reshape(fv.x, [T*C,length(iValid)])';
  Y=reshape(fv.y(2,:), [opt.ncls, nTrials]);
  
  func=['dalmnp300', opt.reg];
  [xx, status]=feval(func, xx0, X, Y, lambda, copy_struct(opt, ...
                                                    'display','eta'), ...
                     'tol', 1e-4, 'solver', 'cg');
  status.dval=status.fval.*(1-status.res);
 otherwise
  if isequal(opt.solver, 'projgrad')
    algo = 'constraint';
  else
    algo = 'squared';
  end


  %% Call the solver
  cdir=cd('opt08');
  fprintf('reg=%s\n',opt.reg);
  problem = packvars(fv.x, ...
                     reshape(fv.y(2,:), [opt.ncls, nTrials]), ...
                     lambda, algo, 'mn', opt.reg);

  [xx, status] = feval(problem.solver, xx0, problem, opt);
  cd(cdir);

end

fprintf('Elapsed time=%g\n', status.time);

cls = struct('W',reshape(xx,[T,C]),'Ws',Ws,'Wt',Wt,'status',status);



