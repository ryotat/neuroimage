% plot_scalps_with_time(mm,nn,tt,W,varargin)
function ax=plot_scalps_with_time(mm,nn,tt,W,mnt,varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'sep', 3,...
                      'clabmrk',[],...
                      'ixt',[]);

ixt = opt.ixt;

if mm==1
  M=1.5;
else
  M = mm+1;
end


if isempty(ixt)
  ixt=sort(pick_indices(sum(abs(W)'),mm*nn,opt.sep));
end

ax.sclp=plot_scalps(M,nn,tt(ixt),W(ixt,:),mnt,'mark_channels',opt.clabmrk);

ax.time=axes('position',[0.3 0.15 0.4 0.9/(mm+2)]);
plot(tt,W);
rr=max(abs(rangeof(W)));
ylim([-1.2*rr 1.2*rr]);
hold on;
plot(ones(2,1)*tt(ixt),ylim'*ones(1,length(ixt)),'b--')
xlabel('Time (ms)');
ytick=get(gca,'ytick');
set(gca,'yticklabel',foreach(@num2str,num2cell(ytick)));



