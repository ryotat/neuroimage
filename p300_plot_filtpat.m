function p300_plot_filtpat(t,C,I,mnt,sgn,varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'fld_W', 'W');

[Uw,Sw,Vw]=svd(getfield(C,opt.fld_W)');
Ws=C.Ws;
Wt=C.Wt;
sd=diag(Sw);

if ~exist('sgn','var') || isempty(sgn)
  sgn=ones(length(sd),1);
end

if size(sgn,1)<size(sgn,2)
  sgn=sgn';
end


if ~iscell(I)
  I=num2cell(I);
end

for ii=1:length(I)
  ix=I{ii}; sdm=mean(sd(ix));
  Sf = Ws*Uw(:,ix)*(sd(ix).*sgn(ix))/sdm;
  Sp = inv(Ws)*Uw(:,ix)*(sd(ix).*sgn(ix))/sdm;
  Tf = Wt*Vw(:,ix)*(sd(ix)./sgn(ix))/sdm;
  Tp = inv(Wt)*Vw(:,ix)*(sd(ix)./sgn(ix))/sdm;

  if length(ix)>1
  stry=sprintf('%d-%d: \\sigma=%s-%s',ix(1),ix(end),num2str_c(sd(ix(1)),2),num2str_c(sd(ix(end)),2));
else
  stry=sprintf('%d: \\sigma=%s',ix,num2str_c(sd(ix),2));
end

  plot_filtpat1(length(I), ii, mnt, Sf, Sp, t, Tp, 'stry', stry);
end

