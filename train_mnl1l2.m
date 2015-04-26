function cls = train_mnl1l2(xTr, yTr, lambda, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'whitening', 'st',...
                      'solver', 'mnl1l2',...
                      'ncls', 6,...
                      'group', 'space',...
                      'W0', [],...
                      'lsmethod', 1,...
                      'display', 'final');

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
nTrials = length(iValid)/opt.ncls;

%% Whitening
[fv,Ws,Wt] = proc_whitenST(fv, opt.whitening);

if strcmp(opt.group,'time')
  fprintf('grouping: time\n');
  xTr = permute(fv.x, [2,1,3]);
elseif strcmp(opt.group,'fourier')
  fprintf('grouping: fourier\n');
  xTr = reshape(fv.x, [1, T*C, length(iValid)]);
  xTr = [real(xTr); imag(xTr)];
else
  fprintf('grouping: space\n');
  xTr = fv.x;
end


%% Call the solver
switch(opt.solver)
  case 'mnl1l2'
   [mm, ix] = max(reshape(fv.y(2,:), [opt.ncls, nTrials]));
   yTr = ix-1;
   time0 = cputime;
   [W, z, status] = feval(opt.solver, xTr, yTr, lambda, opt);
   time1 = cputime;
   status.time = time1-time0;
   fprintf('Elapsed time=%g\n', time1-time0);
 case 'projgrad'
  if size(opt.W0)==[T,C]
    xx0=reshape(opt.W0,[T*C,1]);
  else
    xx0=zeros(T*C,1);
  end
  cdir=cd('opt08');
  problem = packvars(xTr, ...
                   reshape(fv.y(2,:), [opt.ncls, nTrials]), ...
                   lambda, 'constraint', 'mn', 'gl');

  [xx, status] = feval(problem.solver, xx0, problem, opt);
  cd(cdir);
  W=reshape(xx,[T,C])
end

if strcmp(opt.group,'time')
  W=W';
elseif strcmp(opt.group,'fourier')
  Wr = reshape(W(1,:), [T,C]); Wi = reshape(W(2,:), [T,C]);
  W=Wr+sqrt(-1)*Wi;
end
cls = struct('W',W,'Ws',Ws,'Wt',Wt,'status',status);



