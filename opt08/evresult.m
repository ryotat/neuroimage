function evresult(file, ii, varargin)
 
opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'fileres', 'subject_A_result',...
                      'filelog', 'logfile_%g',...
                      'ixn', 1);

file_result = [regexprep(file, '*', opt.fileres), '.mat'];

load(file_result)

if isfield(result, 'C')
  C=result(opt.ixn,ii).C;
else
  C=result(opt.ixn,ii).lambda;
end

file_log    = sprintf([regexprep(file,'*', opt.filelog) '.txt'], C);
S=load(file_log);


subplot(1,2,1);
plot(S(:,1),S(:,5:6),'linewidth',2);
yr=median(S(:,5:6)*[1;-1]);
yc=median(S(:,5));
ylim([yc-yr*5,yc+yr*5]);
grid on;
set(gca,'fontsize',20);
title(untex(file_log));
xlabel('Iterations');
ylabel('Objective value');


subplot(1,2,2);
loglog(S(:,1),abs(S(:,7)),'linewidth',2);
grid on;
set(gca,'fontsize',20);
xlabel('Iterations');
ylabel('Dual penalty');

ax=axes('position',[0.25   0.65   0.17   0.22]);
bar(result(opt.ixn,ii).nm);

set(gcf,'position',[143   468   863   468],...
        'name',[regexprep(file_log(1:end-4),'/','_'),'.eps']);