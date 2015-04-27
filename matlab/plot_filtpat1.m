function ax=plot_filtpat1(mm, ii, mnt, Sf, Sp, tt, Tp, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'mh1', [0.2 0 0.08], ...
                      'mh2', [0.4 0 0.2],...
                      'mv1', [0.1 0 0.2], ...
                      'mv2', [0.2 0.4 0.3],...
                      'amh', [0 0 0],...
                      'amv', [0 0 0],...
                      'N', 4,...
                      'stry', '');

ax(1)=subplotxl(mm, opt.N, 1+opt.N*(ii-1),opt.mh1,opt.mv1,opt);


plotcsp1(mnt, Sf,opt.stry,'fontsize',14);

if ii==1
  title('Filter');
end


ax(2)=subplotxl(mm, opt.N, 2+opt.N*(ii-1),opt.mh1,opt.mv1,opt);
plotcsp1(mnt, Sp);
if ii==1
  title('Pattern');
end

if ~isempty(tt)
  ax(3) = subplotxl(mm, 2, 2+2*(ii-1),opt.mh2,opt.mv2,opt);
  plot(tt, Tp,'linewidth',2);
  pos=get(gca,'position'); pos(3:4) = pos(3:4).*[1.2,1.1];
  % set(gca,'xtick',0:100:600);
  % set(ax,'fontsize',16,'position', pos)
  xlim(rangeof(tt));
  xlabel('Time (ms)');

  grid on;
  if ii==1
    title('Time course');
  end
end