% condor arguments: label, model, lambda, xTrials, file


cdir=pwd;
cd ~/neuro_cvs/matlab/bci
startup_bci
cd(cdir);
% addpath('/home/neuro/mika/CANDY/matlab/condor');
addpath('/home/neuro/ryotat/mutils/');
%addpath('/home/neuro/ryotat/csp/');

 addpath /home/neuro/ryotat/cvx
 addpath /home/neuro/ryotat/cvx/sets
 addpath /home/neuro/ryotat/cvx/keywords
 addpath /home/neuro/ryotat/cvx/builtins
 addpath /home/neuro/ryotat/cvx/commands
 addpath /home/neuro/ryotat/cvx/functions
 addpath /home/neuro/ryotat/cvx/lib
 addpath /home/neuro/ryotat/cvx/structures


% cd /home/neuro/ryotat/SDPT3-4.0-beta/
% startup;
% cd(cdir)


load(file.clab)

%% Load the training data
[cnt, mrk] = p300_preproc(file.dir,...
                          file.train, ...
                          clab, opt.band, opt.fs);

epo = makeEpochs(cnt, mrk, opt.ival);


if opt.bAverage
  fprintf('Averaging subtrials: size(epo.x,3)=%d ->', size(epo.x,3));
  epo = p300_averageSubTrials(epo, opt.nav, opt.nEpochs);
  fprintf('%d\n',size(epo.x,3));
  opt.xvalSubtrials = 1;
else
  epo = p300_sortbycode(epo);
  fprintf('No averaging: size(epo.x,3)=%d\n',size(epo.x,3));
  opt.xvalSubtrials = 15;
end

if ~isempty(xTrials)
  [loss, loss_std]=p300_xvalidation(epo, ...
                                    model, ...
                                    'nEpochs', opt.nEpochs, ...
                                    'nSubtrials', opt.xvalSubtrials,...
                                    'xTrials', xTrials);
else
  loss=[];
  loss_std=[];
end

%% Train the classifier
C = trainClassifier(epo, model.classifier);

%% Load the test data
load(file.true);
[cnt, mrk] = p300_preproc(file.dir,...
                          file.test, ...
                          clab, opt.band, opt.fs, true);

epoTe = makeEpochs(cnt, mrk, opt.ival);
epoTe = p300_sortbycode(epoTe);

%% Apply the classifier
out = applyClassifier(epoTe, model.classifier, C);
out = reshape(out, [opt.nEpochs, opt.nSubtrials, opt.nTrialsTest]);

outcum = cumsum(out,2)./(repmat(1:opt.nSubtrials,[opt.nEpochs,1,opt.nTrialsTest]));

%% Decode and measure the accuracy of ensemble classifier
[fundec, paradec] = getFuncParam(model.decoder);
decodecum = feval(fundec, outcum, paradec{:});
acccum = mean(decodecum==repmat(true,[opt.nSubtrials,1]),2);


vars = {'label',...
        'file',...
        'model',...
        'lambda',...
        'loss',...
        'loss_std',...
        'C',...
        'out*',...
        'decodecum',...
        'acccum',...
        'opt'}

save(file.save, vars{:});


