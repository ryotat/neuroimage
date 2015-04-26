function [result, acccum] = p300_xval(subjects, lambda, model, dir_save, varargin)

global DATA_DIR

opt = propertylist2struct(varargin{:});
opt = set_defaults(opt, 'fs', 60, ...
                        'Nepochs', 12, ...
                        'Nsubtrials', 15, ...
                        'Ntrials', 85, ...
                        'NtrialsTest', 100, ...
                        'band', 20, ...
                        'ival', [0 600], ...
                        'nav', 15, ...
                        'nSets', 1, ...
                        'bAverage', 0,...
                        'basedir', 'P300DATA/',...
                        'dir_data', 'P300DATA/BCI_Comp_III_Wads_2004/',...
                        'label_fmt', 'Subject_%s_lambda_%g',...
                        'file_summ', 1,...
                        'W0', []);


extract(opt);

load([basedir 'clab.mat']);


file = 'Subject_%s_%s';
file_true = [basedir 'test/true_%s.mat'];
file_save = [dir_save label_fmt '.mat'];

if opt.file_summ
  file_summ = [dir_save gen_file_summ(label_fmt,subjects,lambda) '.mat'];
else
  file_summ = '';
end

opt.file_true = file_true;
opt.file_save = file_save;
opt.file_summ = file_summ;


save(file_summ, 'lambda', 'subjects', 'model', 'opt');

for vp=1:length(subjects)
  subject = subjects{vp};

  %% Load the training data
  [cnt, mrk] = p300_preproc(dir_data,...
                            sprintf(file, subject, 'Train'), ...
                            clab, band, fs);
  
  epo = makeEpochs(cnt, mrk, ival);
  
  %% I take the average of 15 repetitions to reduce the number of
  %% training examples from 12*15*85 to 12*85.

  if bAverage
    epo = p300_averageSubTrials(epo, nav, Nepochs);
  else
    epo = p300_sortbycode(epo);
  end

  %% Load the test data
  T=load(sprintf(file_true, subject));
  [cnt, mrk] = p300_preproc(dir_data,...
                            sprintf(file, subject, 'Test'), ...
                            clab, band, fs, T.true);
  
  epoTe = makeEpochs(cnt, mrk, ival);
  epoTe = p300_sortbycode(epoTe);

  save('epoTe.mat', 'epoTe');

  clear cnt mrk epoTe;
  
  model.classifier{3}.W0 = opt.W0;
  for jj=1:length(lambda)
    lmd = lambda(jj);
    fprintf('Subject %s lambda=%g\n\n', subject, lmd);
    model.classifier{2} = lmd;
    C = trainClassifier(epo, model.classifier);
    model.classifier{3}.W0 = C.W;
    
    load('epoTe.mat');
    
    out = reshape(applyClassifier(epoTe, model.classifier, C), [Nepochs, ...
                        Nsubtrials, NtrialsTest, nSets]);
  
    %% cumulative sum of repetitions
    outcumind = cumsum(out,2)./(repmat(1:Nsubtrials,[Nepochs,1,NtrialsTest, nSets]));
  
    %% Decode and measure the accuracy of ensemble classifier
    [fundec, paradec] = getFuncParam(model.decoder);
    decodecum = feval(fundec, mean(outcumind,4), paradec{:});
    acccum = squeeze(mean(decodecum==repmat(T.true,[Nsubtrials,1]),2));

    acccum'

    save(sprintf(file_save,subject,lmd),'C','decodecum','acccum','model','opt');
    
    result(vp,jj) = archive('subject','lmd','C','decodecum','acccum');
    
    clear epoTe;
  end
end

acccum=zeros(Nsubtrials,length(subjects),length(lambda));
for ii=1:prod(size(result))
acccum(:,ii)=result(ii).acccum;
end

if ~isempty(file_summ)
  save(file_summ, 'lambda', 'subjects', 'result', 'acccum', 'model', 'opt');
end



