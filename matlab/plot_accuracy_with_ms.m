function [h,I]=plot_accuracy_with_ms(xvalue, acc, acc_cv, acc_cv_std, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'accmult',100,...
                      'ixrep',[],...
                      'logscale',1,...
                      'ylim','auto',...
                      'ytick','auto',...
                      'fontsize',16,...
                      'mark_lambda',[],...
                      'xlabelstr',[],...
                      'format_ticks',0,...
                      'title','',...
                      'autoselect',0,...
                      'autoselect_n',20,...
                      'colororder',[]);


if opt.logscale
  xx = log(xvalue);
else
  xx = xvalue;
end

if size(acc,1)==1
  acc        = acc';
  acc_cv     = acc_cv';
  acc_cv_std = acc_cv_std';
end


if isempty(opt.ixrep)
  opt.ixrep=1:size(acc,2);
end

if isempty(opt.colororder)
  opt.colororder=get(gca,'colororder');
end

if ~isempty(acc_cv) && all(~isnan(acc_cv(:)))
  ixcv = get_ms_result(acc_cv, opt);
end


if opt.autoselect
  I=findlinspace(xx,opt.autoselect_n);
  I=sort(union(I,ixcv));
else
  I=1:length(xvalue);
end

for kk=1:length(opt.ixrep)
  M=opt.ixrep(kk);
  col = opt.colororder(mod1(kk,size(opt.colororder,1)),:);
  if ~isempty(acc_cv) && all(~isnan(acc_cv(:)))
      % Plot CV curves
      h.cv(kk)=errorbar(xx(I), acc_cv(I,M)*opt.accmult,...
                        acc_cv_std(I,M)*opt.accmult);
      set(h.cv(kk),'color',col,'linewidth',2, 'linestyle', '--');
      hold on;
  end
  h.acc(kk)=plot(xx(I), acc(I,M)*opt.accmult, 'color',col,'linewidth', 2);
  hold on;
end


% Mark MS results
if ~isempty(acc_cv) && all(~isnan(acc_cv(:)))
  plot_ms_result(xvalue, acc, ixcv, opt);
end
axis tight;
if ~isequal(opt.ylim,'auto')
  ylim(opt.ylim);
end

% Mark lambda for visualization
if ~isempty(opt.mark_lambda)
  hold on;
  plot_mark_lambda(xvalue, opt.mark_lambda, opt);
end


set(gca,'fontsize',opt.fontsize)

if opt.logscale
  xrange=roundin(rangeof(xx(abs(xx)~=inf))/log(10));
  xtick = exp((xrange(1):xrange(2))*log(10));
  set(gca,'xtick', log(xtick));
end

if ~isequal(opt.ytick,'auto')
  set(gca,'ytick',opt.ytick);
end

grid on;

if ~isempty(opt.xlabelstr)
  xlabel(opt.xlabelstr);
end

ylabel('Accuracy');

if ~isempty(opt.title)
  title(opt.title);
end

if opt.format_ticks
  format_ticks(gca,num2str_c(xtick));
else
  if opt.logscale
    logticks('x');
  end
end


function ixcv=get_ms_result(acc_cv, opt)

ixcv=zeros(1,length(opt.ixrep));
for kk=1:length(opt.ixrep)
  M=opt.ixrep(kk);
  [mm,ix]=max(acc_cv(:,M));
  ixcv(kk) =medin(find(acc_cv(:,M)==mm));
end

function plot_ms_result(xvalue, acc, ixcv, opt)

colorder=opt.colororder; lenc=size(colorder,1);

if opt.logscale
  xx=log(xvalue);
else
  xx=xvalue;
end

for kk=1:length(opt.ixrep)
  M  = opt.ixrep(kk);
  ix = ixcv(kk);
  if ~isnan(ix)
    plot(xx(ix), acc(ix,M)*opt.accmult, 'o', 'color', ...
         colorder(mod1(kk,lenc),:), 'linewidth', 2, 'markersize', 8);
    fprintf('%s: M=%d ix=%d xval=%g acc=%g\n', opt.title, M, ix, ...
            xvalue(ix), acc(ix,M)*opt.accmult);
  end
  
end

