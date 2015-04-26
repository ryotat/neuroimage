cdir=pwd;
addpath([cdir '/utils']);

pow = -0.25;


subjects = {'A','B'}

Clist=1./exp(linspace(log(1),log(1e+4),20));


Tab = cell2mat({'A','B','C','D','E','F';...
       'G','H','I','J','K','L';...
       'M','N','O','P','Q','R';...
       'S','T','U','V','W','X';...
       'Y','Z','1','2','3','4';...
       '5','6','7','8','9','_'});


clsopt = struct('W0',[],'nSamples', 37, 'nChannels', 64, 'ncls', 6, ...
                'group','space','solver','projgrad', 'whitening', 'scalest', 'display', 2);
 
classy = {'mnl1l2', [], clsopt};
           

model = struct('classifier', {classy},...
              'decoder', {{'p300_decode', Tab}});

dir_save = 'P300DATA/results/space/';

label_fmt = sprintf('Subject=%%s_C=%%g_space',pow);

[result,acccum]=p300_xval(subjects,...
                          Clist,...
                          model,...
                          dir_save,...
                          'label_fmt', label_fmt);

