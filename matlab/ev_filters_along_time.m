% lambda=exp(linspace(log(0.01),log(100),20));
% dir_result = 'P300DATA/results/mnl1l2_oldb/';
% label_fmt='Subject_%s_lambda_%g_space_scalest'
% load mnt_p300.mat
% t=0:1000/60:600;
% slideShow('ev_filters_along_time',[1 20],'dir_result',dir_result,'label_fmt',label_fmt,'lambda',lambda,'group','','ixt',[],'t',t,'mnt',mnt,'ev_title','','sbj','B','rows',1,'sep',1);

if ~exist('rows','var')
  rows=1;
end

if ~exist('cols','var')
  cols=7;
end

if ~exist('sep','var')
  sep=3;
end

if ~exist('clabmrk','var')
  clabmrk=[];
end

if ~exist('pat','var')
  pat = [];
end


if isempty(group)
  file=[dir_result sprintf(label_fmt, sbj, lambda(i)), '.mat'];
  
  if ~isempty(findstr(label_fmt, 'space'))
    group='space';
  elseif ~isempty(findstr(label_fmt, 'time'))
    group='time';
  else
    group='ds';
  end
else  
  file=[dir_result sprintf(label_fmt, sbj, lambda(i), group), '.mat'];
end

S=load(file);

if isfield(S,'acccum')
  acc = S.acccum([5 15]);
else
  acc = S.acc;
end


C=S.C;
if ~isfield(C,'W')
  C.W=C.W1;
end


if isempty(pat) || pat==0
  W=real(C.Wt*C.W*C.Ws);
else
  W=real(C.Wt^(-pat(1))*C.W*C.Ws^(-pat(2)));
end

ax=plot_scalps_with_time(rows,cols,t,W,mnt,'ixt',ixt,'sep',sep,'mark_channels',clabmrk);


if rows==1
  M=1.5;
else
  M = rows+1;
end

set(gcf,'position',[4 30 1024 232*M]);
  
set(gcf,'papersize',[20 11]);


if length(acc)>1
str=sprintf('sbj=%s lmd=%g nm=%g (%s) acc=[%g %g]',sbj,lambda(i), ...
            normc(C.W,group), group, round(acc([1 end])'*100))
else
str=sprintf('sbj=%s lmd=%g nm=%g (%s) acc=[%g]',sbj,lambda(i), ...
            normc(C.W,group), group, acc)
end
fprintf('%s_%g_%s.eps\n',sbj,lambda(i),group);

set(gcf,'name',str);

axes(ax.sclp(1));
text(0,1.3,[file '  ' str],'unit','normalized','interpreter','none');
