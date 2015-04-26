function cls = train_ds_1st2nd(xTr, yTr, lambda, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'whitening', 'st',...
                      'loss','lr',...
                      'reg', 'ds',...
                      'solver', 'projgrad',...
                      'ncls', 2,...
                      'W0', [],...
                      'display', 1,...
                      'epsf',1e-4,...
                      'maxiter',100000,...
                      'adapt_lambda', 1,...
                      'eta',1);

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
yTr=yTr(:,iValid);

if ndims(xTr)<3
  T = opt.nSamples;
  C = opt.nChannels;
  n = length(iValid);

  if sz(1)==T*C
    fv   = struct('x',reshape(xTr, [T, C, n]),'y',yTr);
    fvcov= [];
  else
    [X1,X2]=unpack1st2nd(xTr(:,iValid),T,C,n,'flat',1);
    fv   = struct('x',X1,'y',yTr);
    fvcov= cell(1,length(X2));
    for ii=1:length(X2)
      fvcov{ii}=struct('x',X2{ii},'y',yTr);
    end
  end
else
  [T,C,n]=size(xTr(:,:,iValid));
  fv    = struct('x', xTr(:,:,iValid), 'y', yTr(:,iValid));
  fvcov = {proc_covariance(fv)};
end


if size(yTr,1)~=opt.ncls
  nTrials = length(iValid)/opt.ncls;
else
  nTrials = length(iValid);
end

%% Whitening
[fv,Ws,Wt] = proc_whitenST(fv, opt.whitening);

if ~isequal(opt.whitening, 'none')
  Ws2=cell(1,length(fvcov));
  for ii=1:length(fvcov)
    [fvcov{ii},Ws2{ii}]=proc_whiten(fvcov{ii});
  end
else
  Ws2=repmat({eye(C)},[1,length(fvcov)]);
end

X = pack1st2nd(fv, fvcov);


switch (opt.solver)
 case 'projgrad'
  algo = 'constraint';
  if opt.adapt_lambda>0
    opt.adapt_lambda=1;
  end
 case 'sublbfgs'
  algo = 'squared';
  if opt.adapt_lambda>0
    opt.adapt_lambda=2;
  end
 otherwise
  %error('Unknown solver [%s].',opt.solver);
end

%% Initial value
if ~isempty(opt.W0)
  if prod(size(opt.W0))==size(X,1)
    xx0=reshape(opt.W0,[size(X,1),1]);
  else
    xx0=opt.W0;
  end
else
  xx0=zeros(size(X,1),1);
end

%% Call the solver
if length(lambda)==1 && opt.adapt_lambda>0
  ftmp=inline('sum(var(reshape(fv.x,[size(fv.x,1)*size(fv.x,2),size(fv.x,3)])''));','fv');
  vv=cell2mat(foreach(ftmp,[{fv},fvcov]));
end


if isequal(opt.solver, 'dal')
  Y=(2*fv.y(2,:)-1)';
  fv=[{fv} fvcov];
  [X,ns,nc]=pack_matrices(fv,'weight',sqrt(1./vv));
  X=X'; blks = [ns; nc]';
  [ww,bias,status]=dallrds(zeros(size(X,2),1),1e-6*randn(1),X,Y,lambda,'blks',blks,'eta',opt.eta,'tol',opt.epsf,'display',opt.display); 

  status.dval=status.fval.*(1-status.res);
  
  xx=[ww;bias];

  [W1,W2]=unpack1st2nd(ww,T,C,1,'weight',sqrt(1./vv));
  cls = strukt('W1',W1,'W2',W2,'bias',bias,'Ws',Ws,'Wt',Wt,'Ws2',Ws2,'status',status);

else
cdir=cd('opt08');
  fprintf('loss=%s reg=%s\n',opt.loss, opt.reg);
  if size(fv.y,1)~=opt.ncls
    Y=reshape(fv.y(2,:), [opt.ncls, nTrials]);
  else
    Y=2*fv.y(2,:)-1;
  end
  
  %% determine lambda
  if length(lambda)==1 && opt.adapt_lambda>0
    lambda = lambda*((vv/vv(1)).^(opt.adapt_lambda/2))
  end
  
  
  problem = packvars([{fv},fvcov], Y, lambda, algo, opt.loss, opt.reg);

  if problem.hasbias && length(xx0)<size(problem.X,1)
    xx0 = [xx0;0];
  end
  
  [xx, status] = feval(problem.solver, xx0, problem, opt);

  cd(cdir);

  if problem.hasbias
  [W1,W2]=unpack1st2nd(xx(1:end-1),T,C,1);
  cls = strukt('W1',W1,'W2',W2,'bias',xx(end),'Ws',Ws,'Wt',Wt,'Ws2',Ws2,'status',status);
  else
  [W1,W2]=unpack1st2nd(xx,T,C,1);
  cls = strukt('W1',W1,'W2',W2,'Ws',Ws,'Wt',Wt,'Ws2',Ws2,'status',status);
  end

end


cls.v = vv;

fprintf('Elapsed time=%g\n', status.time(end));





