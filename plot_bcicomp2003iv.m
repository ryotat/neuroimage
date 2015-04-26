function [h,acc,D,I]=plot_bcicomp2003iv(lambda, memo, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'reg','ds',...
                      'normx',0,...
                      'v', [],...
                      'div',0,...
                      'fld_C','C',...
                      'fld_acc_cv','loss',...
                      'fld_acc_cv_std','loss_std',...
                      'fld_acc','acc',...
                      'color','b',...
                      'autoselect',0,...
                      'autoselect_n',20);

acc_cv=[];
acc_cv_std=[];
acc=[];

if ~isempty(opt.fld_acc_cv)
  try
    acc_cv=cell2mat(getfieldarray(memo,opt.fld_acc_cv));
  end
end

if ~isempty(opt.fld_acc_cv_std)
  try
    acc_cv_std=cell2mat(getfieldarray(memo,opt.fld_acc_cv_std));
  end
end

if ~isempty(opt.fld_acc)
  try
    acc=cell2mat(getfieldarray(memo,opt.fld_acc));
  end
end


% nmtmp = normc(getfield(memo(1),opt.fld_C),opt.reg,opt,'dosum',0);
% D=zeros(length(memo),length(nmtmp));
for ii=1:length(memo)
  C=getfield(memo(ii),opt.fld_C);
  D(ii,:)=normc(C,opt.reg,opt,'dosum',0);
  switch(opt.reg)
   case 'ds'
    nm(ii)=sum(D(ii,:));
   case 'fro'
    nm(ii)=sqrt(sum(D(ii,:).^2));
  end
end    


if opt.normx
  xvalue = nm;
  switch(opt.reg)
   case 'ds'
     xlabelstr = '\Omega_{DS}(\theta)';
   case 'fro'
    xlabelstr = '\Omega_F(\theta)';
  end
else
  xvalue = lambda;
  xlabelstr = 'Regularization constant \lambda';
end


[h,I]=plot_accuracy_with_ms(xvalue, acc, acc_cv, acc_cv_std,...
                      opt, ...
                      'lambda', lambda, ...
                      'xlabelstr', xlabelstr, ...
                      'colororder', opt.color,...
                      'accmult', 1);

% $$$ mm=max(acc_cv);ixcv=medin(find(acc_cv==mm));
% $$$ fprintf('lambda=%g acc_cv=%g acc=%g\n', lambda(ixcv), acc_cv(ixcv), acc(ixcv));
% $$$ 
% $$$ 
% $$$ if opt.autoselect
% $$$   I=findlinspace(log(nm),opt.autoselect_n);
% $$$   I=sort(union(I,ixcv));
% $$$ else
% $$$   I=1:length(xvalue);
% $$$ end
% $$$ 
% $$$ 
% $$$ if ~isempty(acc_cv)
% $$$   h.cv=errorbar(log(xvalue(I)),acc_cv(I),acc_cv_std(I));
% $$$   set(h.cv,'color',opt.col_cv,'linestyle',opt.style_cv,'linewidth',2);
% $$$ end
% $$$ hold on;
% $$$ 
% $$$ h.acc=plot(log(xvalue(I)),acc(I),'color',opt.col_acc,'linestyle',opt.style_acc,'linewidth',2);
% $$$ 
% $$$ grid on;
% $$$ 
% $$$ plot(log(xvalue(ixcv)), acc(ixcv), 'o', 'color',opt.col_acc,'linewidth', 2, 'MarkerSize',8);
% $$$ 
% $$$ logticks('x');
% $$$ xlabel(xlabelstr);
% $$$ ylabel('Accuracy');
% $$$ axis tight;
% $$$ set(gca,'ylimmode','auto');


