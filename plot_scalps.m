function ax=plot_scalps(mm,nn,tt,W,mnt,varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt,'mark_channels',[]);

r=max(real(rangeof(W)));
for jj=1:length(tt)
  ax(jj)=subplotxl(mm,nn,jj, [0.01 0.01 0.05], 0.01);
  scalpPlot(mnt, W(jj,:),'colAx',[-.5*r .5*r],'scalePos','none','mark_channels',opt.mark_channels);
  title(sprintf('t=%d',round(tt(jj))),'fontsize',16);
end
