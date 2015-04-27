function [memo,opt] = xval_bcicomp2003iv_mb(lambda, model, dir_result, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'band1', 20, 'band2', {},...
                      'file_suf', '',...
                      'datadir','data/bcicomp2003iv/',...
                      'fs_ogn', 1000,...
                      'fs', 100,...
                      'xTrials', []);

S=load([opt.datadir 'sp1s_aa_1000Hz.mat']);
load([opt.datadir 'labels_data_set_iv.txt']);

y = [1-S.y_train; S.y_train];
yte = [1-labels_data_set_iv'; labels_data_set_iv'];


fprintf('Extracting 1st order component band=[%d]\n', opt.band1);
x1_tr = filtButterDownSample(S.x_train, 5, opt.band1, opt.fs_ogn, opt.fs_ogn/opt.fs);
x1_te = filtButterDownSample(S.x_test, 5, opt.band1, opt.fs_ogn, opt.fs_ogn/opt.fs);
fv1_tr = struct('x', x1_tr, 'y', y,...
                'clab', {S.clab}, 'fs', opt.fs, 't', -620:10:-130,...
                'className', {{'left','right'}});
fv1_te = struct('x', x1_te, 'y', yte,...
                'clab', {S.clab}, 'fs', opt.fs, 't', -620:10:-130,...
                'className', {{'left','right'}});

if ~isempty(opt.band2)
  fv2_tr=cell(1,length(opt.band2));
  fv2_te=cell(1,length(opt.band2));
  for ii=1:length(opt.band2)
    fprintf('Extracting 2nd order component band=[%d]\n', opt.band2{ii});
    x2_tr = filtButterDownSample(S.x_train, 5, opt.band2{ii}, opt.fs_ogn, opt.fs_ogn/opt.fs);
    x2_te = filtButterDownSample(S.x_test, 5, opt.band2{ii}, opt.fs_ogn, opt.fs_ogn/opt.fs);
    fv2_tr{ii} = proc_covariance(struct('x', x2_tr, 'y', y,...
                'clab', {S.clab}, 'fs', opt.fs, 't', -620:10:-130,...
                'className', {{'left','right'}}));
    fv2_te{ii} = proc_covariance(struct('x', x2_te, 'y', yte,...
                'clab', {S.clab}, 'fs', opt.fs, 't', -620:10:-130,...
                'className', {{'left','right'}}));
  end
else
  fv2_tr={};
  fv2_te={};
end

if length(opt.band2)>0
  fvTr = proc_pack1st2nd(fv1_tr, fv2_tr{:});
  fvTe = proc_pack1st2nd(fv1_te, fv2_te{:});
else
  fvTr = fv1_tr;
  fvTe = fv1_te;
end

keyboard;

xx0=model{3}.W0;
for ii=1:length(lambda)
  lmd = lambda(ii);
  model{2} = lmd;

  model{3}.W0 = xx0;
  C   = trainClassifier(fvTr, model);
  xx0 = vectc(C);
  out = applyClassifier(fvTe, model, C);
  acc =  100*(1-mean(loss_0_1(fvTe.y, out)));

  
  % model selection
  model{3}.W0 = [];
  if ~isempty(opt.xTrials)
    [loss, loss_std]=xvalidation(fvTr, model, 'xTrials', opt.xTrials);
    loss = (1-loss)*100;
    loss_std = loss_std*100;
  else
    loss=[]; loss_std=[];
  end
  
  nm = normc(C,'ds','v',C.v,'div',0);

  
  fprintf('lambda=%g (nm=%g) fval=%g dval=%g acc=%g\n', lambda(ii), ...
          nm, C.status.fval(end), C.status.dval(end), acc);

  memo(ii)=archive('lmd','model','C','out','loss','loss_std','acc','nm');
  saveStruct([dir_result sprintf('lambda=%g',lmd) opt.file_suf '.mat'],memo(ii));
end

