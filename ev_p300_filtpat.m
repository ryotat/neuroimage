% slideShow('ev_p300_filtpat',[1 20],'lambda',Clist,'t',t,'mnt',mnt,'dir_result',dir_result,'label_fmt',label_fmt,'I',1:2,'sbj','A','ev_title','')

if ~exist('sgn','var')
  sgn=[];
end


S=load([dir_result sprintf(label_fmt, sbj, lambda(i)) '.mat']);
str=sprintf('[%d] sbj=%s lmd=%g nm=%g acc=[%g %g]\n',i,sbj,lambda(i), normc(S.C.W,'ds'), ...
        S.acccum(5), S.acccum(15));


p300_plot_filtpat(t,S.C,I,mnt,sgn);

if length(I)==2
  set(gcf,'position',[69   154  663   358]);
else
  set(gcf,'position',[69   154   595   512]);
end


fprintf(str);
set(gcf,'name',str);
text(1.1,0,str,'unit','normalized','rotation',90,'interpreter','none');
