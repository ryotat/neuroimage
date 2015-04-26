% p300_xvalidation - xvalidates P300 data with respect to the
%                    character prediction accuracy
%
% Syntax:
%  [loss, loss_std, out_test, memo] = p300_xvalidation(epo, model, varargin)
%
%
% Ryota Tomioka 2009

function [loss, loss_std, out_test, memo] = p300_xvalidation(epo, model, varargin)
  
opt = propertylist2struct(varargin{:});
opt = set_defaults(opt, 'nEpochs', 12,...
                        'nSubtrials', 1,...
                        'xTrials', [10, 10],...
                        'sample_fcn', 'kfold',...
                        'divTr', [],...
                        'divTe', [],...
                        'dispSubtrials', [5 15],...
                        'averageloss', 1);

opt.dispSubtrials = unique(min(opt.dispSubtrials,opt.nSubtrials));


nTrials = size(epo.y,2)/opt.nEpochs/opt.nSubtrials;

if isempty(opt.divTr) || isempty(opt.divTe)
  sample_fcn = ['sample_' opt.sample_fcn];
  [divTr, divTe] = feval(sample_fcn, ones(1,nTrials), opt.xTrials);
else
  divTr = opt.divTr;
  divTe = opt.divTe;
end

[fundec, paradec] = getFuncParam(model.decoder);


out_test = zeros([length(divTe), opt.nEpochs*opt.nSubtrials*nTrials]);

loss = zeros(opt.nSubtrials, length(divTe),length(divTe{1}));

for ii=1:length(divTe)
  for jj=1:length(divTe{ii})
    nTrain = length(divTr{ii}{jj});
    nTest  = length(divTe{ii}{jj});
    
    fprintf('Fold (%d,%d): nTrain=%d nTest=%d\n',ii,jj,nTrain,nTest);
    [epoTr, indexTr] = p300_selectTrials(epo, divTr{ii}{jj},...
                              opt.nEpochs, opt.nSubtrials, nTrials);
    
    [epoTe, indexTe] = p300_selectTrials(epo, divTe{ii}{jj},...
                              opt.nEpochs, opt.nSubtrials, nTrials);
    
    C = trainClassifier(epoTr, model.classifier);
    out = applyClassifier(epoTe, model.classifier, C);
    outcum=cumsum(reshape(out, [opt.nEpochs, opt.nSubtrials, nTest]),2)./(repmat(1:opt.nSubtrials,[opt.nEpochs,1,nTest]));

    
    yav   = squeeze(mean(reshape(epoTe.y(2,:), [opt.nEpochs,opt.nSubtrials, nTest]),2));
    
    decodecum = feval(fundec, outcum, paradec{:});
    true   = feval(fundec, yav, paradec{:});
    
    loss(:, ii,jj)=1-mean(decodecum==repmat(true,[opt.nSubtrials,1]),2);
    
    out_test(ii, indexTe) = out;

    memo(ii,jj) = archive('indexTr', 'indexTe', 'C', 'out', 'decodecum','true');
    
    if opt.nSubtrials==1
      fprintf('true   =%s\ndecoded=%s\naccuracy=%g%%\n', true, decode, (1-loss(ii,jj))*100);
    else
      for kk=1:length(opt.dispSubtrials)
        fprintf('[%02d] decoded=%s : %g%%\n', opt.dispSubtrials(kk), decodecum(opt.dispSubtrials(kk),:), 100*(1-loss(opt.dispSubtrials(kk),ii,jj)));
      end
      fprintf('---------------------------------------\n     true   =%s\n',true);
    end
  end
end


loss_std = std(mean(loss,3),[],2);

if opt.averageloss
  loss = mean(mean(loss,3),2);
end



    

