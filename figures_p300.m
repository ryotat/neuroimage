dir = './';

t=0:1000/60:600;

load([dir 'mnt_p300.mat'])

dir_result=[dir 'P300DATA/results/mnds_fro/'];
label_fmt='Subject=%s_lambda=%g_exp=-0.25'
lambda=exp(linspace(log(1e+6),log(1),20));
plot_result_p300(dir_result, label_fmt, 'showspec',1,'lambda',lambda,'ylim',[40 100],'normx',1,'reg','fro');


dir_result=[dir 'P300DATA/results/mnds/'];
label_fmt='Subject=%s_C=%g_exp=-0.25';
Clist=1./exp(linspace(log(0.1),log(10),20));
lmdmark=[1.4385    1.1288];
[D,memo]=plot_result_p300(dir_result, label_fmt, 'showspec',1,'lambda',Clist,'ylim',[40 100],'normx',1,'mark_lambda',lmdmark);

IA=get_ms_result(-memo(1).loss, 1:15);
IB=get_ms_result(-memo(2).loss, 1:15);

S=load([dir 'P300DATA/raktom/result4.mat']);

%M =[1, 2, 3, 4, 5, 10,13,15];
%RG=[16,32,52,60,72,83,94,97;...
%    35,53,62,68,75,91,96,96];


figure     

subplot(2,1,1);
plot(1:15, memo(1).acccum(IA+(0:14)*20)*100,'-o', 1:15, S.perfA,'-x', ...
      'linewidth',2)
set(gca,'fontsize',16);
grid on;
ylabel('Accuracy');
title('Subject A');
h=legend('DS regularization','R&G');
set(h,'fontsize',12);
subplot(2,1,2);
plot(1:15, memo(2).acccum(IB+(0:14)*20)*100,'-o', 1:15, S.perfB,'-x', ...
      'linewidth',2)
set(gca,'fontsize',16);
grid on;
xlabel('Number of repetitions M');
ylabel('Accuracy');
title('Subject B');
h=legend('DS regularization','R&G');
set(h,'fontsize',12);


% subject A

% subject B
slideShow('ev_p300_filtpat',10,'lambda',Clist,'t',t,'mnt',mnt,'dir_result',dir_result,'label_fmt',label_fmt,'I',1:3,'sbj','B','ev_title','','sgn',[-1 1 1])

slideShow('ev_p300_filtpat',10,'lambda',Clist,'t',t,'mnt',mnt,'dir_result',dir_result,'label_fmt',label_fmt,'I',{1:2},'sbj','B','ev_title','','sgn',[-1 1 1])

% Channel selection regularizer
lambda=exp(linspace(log(0.01),log(100),20));
dir_result = [dir 'P300DATA/results/mnl1l2_oldb/'];
label_fmt='Subject_%s_lambda_%g_space_scalest'
lmdmark=[8.86 8.86];
[D,memo]=plot_result_p300(dir_result, label_fmt,...
                 'showspec',1,'normx',1,'reg','gl','group','space','ylim',[40 100],'mark_lambda',lmdmark);

% Subject A
slideShow('ev_filters_along_time',15,'dir_result',dir_result,'label_fmt',label_fmt,'lambda',lambda,'group','','ixt',[],'t',t,'mnt',mnt,'ev_title','','sbj','A','rows',1,'sep',1);

% Subject B
 slideShow('ev_filters_along_time',15,'dir_result',dir_result,'label_fmt',label_fmt,'lambda',lambda,'group','','ixt',[],'t',t,'mnt',mnt,'ev_title','','sbj','B','rows',1,'sep',2);

%% Temporal basis selection regularizer
label_fmt='Subject_%s_lambda_%g_time_scalest'
lmdmark = [8.86 14.4];
plot_result_p300(dir_result, label_fmt,...
 'showspec',1,'normx',1,'reg','gl','group','time','ylim',[40 100],'mark_lambda',lmdmark,'ytick',40:10:100);

% Subject A
slideShow('ev_filters_along_time',15,'dir_result',dir_result,'label_fmt',label_fmt,'lambda',lambda,'group','','ixt',[],'t',t,'mnt',mnt,'ev_title','','sbj','A','rows',1,'sep',1);

% Subject B
slideShow('ev_filters_along_time',16,'dir_result',dir_result,'label_fmt',label_fmt,'lambda',lambda,'group','','ixt',[],'t',t,'mnt',mnt,'ev_title','','sbj','B','rows',1,'sep',1);

%label_fmt='Subject_%s_lambda_%g_fourier_scales+fourier'
% plot_result_p300(dir_result, label_fmt,...
% 'showspec',1,'normx',1,'reg','gl','group','fourier','ylim',[40 100]);
