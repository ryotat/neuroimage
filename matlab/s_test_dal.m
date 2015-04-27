load P300DATA/clab.mat

band = 20;
fs   = 60;
ival = [0 600];

pow = -0.25;


subjects = {'A'}; %{'A','B'}

Clist=1/5; % 1./exp(linspace(log(0.1),log(10),20));


Tab = cell2mat({'A','B','C','D','E','F';...
       'G','H','I','J','K','L';...
       'M','N','O','P','Q','R';...
       'S','T','U','V','W','X';...
       'Y','Z','1','2','3','4';...
       '5','6','7','8','9','_'});


clsopt = struct('W0',[],'nSamples', 37, 'nChannels', 64, 'ncls', 6, ...
                'solver','dal', 'whitening', {{'st',pow,pow}}, 'display', 3);
 
classy = {'mnds', [], clsopt};
           

model = struct('classifier', {classy},...
              'decoder', {{'p300_decode', Tab}});


dir_data='P300DATA/BCI_Comp_III_Wads_2004/';
file = 'Subject_%s_%s';
subject='A';

[cnt, mrk] = p300_preproc(dir_data,...
                          sprintf(file, subject, 'Train'), ...
                          clab, band, fs);

epo = makeEpochs(cnt, mrk, ival);
epo = p300_sortbycode(epo);

classy{2}=10;
C = trainClassifier(epo, classy);

