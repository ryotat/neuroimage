function [D,memo]=plot_result_p300(dir_result, label_fmt, varargin)

global DATA_DIR;

opt = propertylist2struct(varargin{:});
opt = set_defaults(opt, 'reg','ds',...
                        'group','space',...
                        'ixrep',[5,15],...
                        'ixfold', 4,...
                        'bShowLikelihood', 0,...
                        'subjects', {'A','B'},...
                        'lambda', exp(linspace(log(0.01), log(100), 20)),...
                        'showspec',1,...
                        'normx',0,...
                        'thzero', 0.01,...
                        'plotbar', 0,...
                        'nSamples', 37,...
                        'nChannels', 64,...
                        'ylim', [50 100],...
                        'ytick', 50:10:100,...
                        'fontsize',16,...
                        'mark_lambda', [nan nan]);



extract(opt);


% dir_data = [DATA_DIR 'eegImport/bci_competition_iii/albany/'];
file = 'Subject_%s_%s';

load('P300DATA/clab.mat');

fs = 60;
band = 20;
ival  = [0 600];
nav = 15;

Nepochs    = 12;
Nsubtrials = 15;
Ntrials    = 85;
NtrialsTest= 100;


Tab = cell2mat({'A','B','C','D','E','F';...
       'G','H','I','J','K','L';...
       'M','N','O','P','Q','R';...
       'S','T','U','V','W','X';...
       'Y','Z','1','2','3','4';...
       '5','6','7','8','9','_'});

switch(opt.reg)
 case 'ds',
   D=zeros(length(lambda), min(nSamples,nChannels), length(subjects));
 case 'gl',
  switch(opt.group)
   case 'time',
    D=zeros(length(lambda), nSamples, length(subjects));
   case 'space',
    D=zeros(length(lambda), nChannels, length(subjects));
   case 'fourier',
    D=zeros(length(lambda), nSamples*nChannels, length(subjects));
  end
end


nRows = 1+(opt.showspec>0);

figure;
for ii=1:length(subjects)
  subject = subjects{ii};
  
  [loss, loss_std, acccum, D(:,:,ii)]=load_result_p300(dir_result, label_fmt, subject, opt);
  
  memo(ii)=archive('subject','loss','loss_std','acccum');
  
  if opt.normx
    if isequal(opt.reg,'fro')
      xvalue = sqrt(sum(D(:,:,ii).^2,2))';
      xlabelstr = '    \Omega_F(\theta)';
    else
      xvalue = sum(D(:,:,ii),2)';
      switch(opt.reg)
       case 'ds'
         xlabelstr = '    \Omega_{DS}(\theta)';
       case 'gl'
        switch(opt.group)
         case 'space'
          xlabelstr = '    \Omega_{C}(\theta)';
         case 'time'
          xlabelstr = '    \Omega_{T}(\theta)';
         case 'fourier'
          xlabelstr = '    \Omega_{SSC}(\theta)';
        end
      end
    end
  else
    xvalue = lambda;
    xlabelstr = 'Regularization constant \lambda';
  end

  
  % Plot performance
  subplotxl(nRows,2,ii,[0.12 0.12 0.05], [0.15 0.12 0.1]);
  
  plot_accuracy_with_ms(xvalue, acccum, ...
                        1-loss, loss_std,...
                        opt,...
                        'xlabelstr',xlabelstr,...
                        'title',['Subject ' subject],...
                        'mark_lambda',opt.mark_lambda(ii),...
                        'format_ticks',1,...
                        'logscale',1);

  xlim1  = get(gca,'xlim');
  xtick1 = get(gca,'xtick');
  
  % Plot counts/spectrum
  if opt.showspec>0
    subplotxl(nRows, 2,ii+2,[0.12 0.12 0.05],[0.15 0.12 0.1]);
    N=countComp(D(:,:,ii)','tolrel',opt.thzero)';
    switch(opt.showspec)
     case 1
      plot(log(xvalue), N,'linewidth',2);
      ylim([0, ceil(max(N)/10)*10]);

      if isequal(opt.reg,'gl') && isequal(opt.group,'space')
        set(gca,'ytick',[0 16 32 48 64]);
      else
        set(gca,'ytick',[0 10 20 30 40]);
      end
      
      ylabel('# active components')
     case 2
      plot(log(xvalue), D(:,:,ii),'linewidth',2);
      ylim([0, max(median(D(:,:,ii)))*2]);
      switch(opt.reg)
       case {'ds','fro'}
        ylabel('Singularvalues');
       case'gl'
        ylabel('Component norms');
      end
    end      
    if ~isempty(opt.mark_lambda)
      hold on;
      plot_mark_lambda(xvalue, opt.mark_lambda(ii), opt);
    end
    

    set(gca,'fontsize', opt.fontsize,...
            'xtick'   , xtick1,...
            'xlim'    , xlim1);
    grid on;
    xlabel(xlabelstr);
    format_ticks(gca,num2str_c(exp(xtick1)));
  end
  
end

text(1.1,0,[dir_result label_fmt],'unit','normalized','rotation',90,'interpreter','none');

if opt.showspec>0
  set(gcf,'position',[164   326   600   546]);
else
  set(gcf,'position',[164   332   713   334]);
end

if opt.plotbar
  slideShow('bar(D(i,:,1))',[1,20],'D',D);
end




function [loss, loss_std, acccum, D]=load_result_p300(dir_result, label_fmt, subject, opt)

lambda = opt.lambda;

switch(opt.reg)
 case {'ds' 'fro'}
   D=zeros(length(lambda), min(opt.nSamples,opt.nChannels));
 case 'gl',
  switch(opt.group)
   case 'time',
    D=zeros(length(lambda), opt.nSamples);
   case 'space',
    D=zeros(length(lambda), opt.nChannels);
   case 'fourier',
    D=zeros(length(lambda), opt.nSamples*opt.nChannels);
  end
end

loss=nan*ones(length(lambda), 15);
loss_std=nan*ones(length(lambda), 15);
acccum=nan*ones(length(lambda), 15);

for jj=1:length(lambda)
  label = sprintf(label_fmt, subject, lambda(jj));
  file = [dir_result label '.mat'];
  try
    S=load(file);
  catch
    fprintf('file[%s] not found...\n',file);
    continue;
  end
  
  switch(S.model.classifier{1})
   case 'p300_div'
    C = S.C.cls(opt.ixfold);
   case {'mnl1l2','mnds'}
    C = S.C;
  end
  
  if isfield(S,'loss') && ~isempty(S.loss)
    loss(jj,:)=S.loss';
    loss_std(jj,:)=S.loss_std';
  else
    % loss=[];
    % loss_std=[];
  end
  acccum(jj,:) = S.acccum';

  switch(opt.reg)
   case {'ds','fro'}
    D(jj,:) = normc(C, opt.reg, 'dosum', 0);
   case 'gl'
    D(jj,:) = normc(C, opt.group, 'dosum', 0);
  end
  
end % end loop over lambda

function ix=plot_mark_lambda(xvalue, lmd, opt)

ix = findnearest(opt.lambda, lmd);
if ~isempty(ix)
  plot(ones(2,1)*log(xvalue(ix)), ylim'*ones(1,length(ix)),'--', 'color',[.5 .5 .5],'linewidth',2);
end

fprintf('lmd=%g ix=%d nm=%g\n', lmd, ix, xvalue(ix));
