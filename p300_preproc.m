% p300_preproc - Preprocess P300 dataset from BCI competition III
% 
% Syntax:
%  [cnt,mrk] = p300_preproc(datadir, dataset, clab, band, fs, TargetChar)
%
% Ryota Tomioka 2009

function [cnt_flt, mrk] = p300_preproc(datadir, dataset, clab, band, fs, TargetChar)

  Tab = cell2mat({'A','B','C','D','E','F';...
       'G','H','I','J','K','L';...
       'M','N','O','P','Q','R';...
       'S','T','U','V','W','X';...
       'Y','Z','1','2','3','4';...
       '5','6','7','8','9','_'});

  
  file = [datadir dataset];
  load(file);
  [Ntrials, Ntmp, Nchannels] = size(Signal);
  
  Signal   = double(reshape(permute(Signal, [2,1,3]), [Ntmp*Ntrials, ...
                      Nchannels]));
  
  Flashing = double(reshape(Flashing', [1,Ntmp*Ntrials]));
  StimulusCode = double(reshape(StimulusCode', [1,Ntmp*Ntrials]));
  
  pos = find([1, diff(Flashing)]==1);
  code = StimulusCode(pos);

  if exist('TargetChar','var')
    y = zeros(1, length(pos));
    for ii=1:length(TargetChar)
      ix = (1:180) + (ii-1)*180;
      [irow, icol] = ind2sub(size(Tab), find(Tab==TargetChar(ii)));
      y(ix) = double(code(ix)==icol | code(ix)==irow+6);
    end
    
    %% Check
    codey=code(find(y));
    clear L;
    for ii=1:length(TargetChar)
      cor = unique(codey((1:30)+(ii-1)*30));
      L(ii)=Tab(cor(2)-6, cor(1));
    end
    if ~isequal(L, TargetChar)
      fprintf('y is corrupted\n')
      keyboard;
    end
    
    y = [1-y; y];
  else
    y = [];
  end
    
  cnt = strukt('x', Signal, ...
               'fs', 240, ...
               'clab', clab, ...
               'file', [file '.mat'], ...
               'title', ['Data Set II - ' dataset]);
  
  mrk = strukt('pos', pos, ...
               'fs', 240, ...
               'code', code,...
               'y', y,...
               'indexedByEpochs', {'code'},...
               'className', {'standard', 'deviant'});

  [b,a]=getbutter(5, band, cnt.fs);
  cnt_flt = proc_filt(cnt, b, a);

  [cnt_flt, mrk] = proc_downsample(cnt_flt, mrk, ceil(cnt.fs/fs));
