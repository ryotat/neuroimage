function [Xout,dd1,dd2,n]=pack_matrices(X, varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'weight',ones(1,length(X)));

if length(opt.weight)~=length(X)
  error('Opt.weight must have the same length as X.');
end

dd1=zeros(1,length(X));
dd2=zeros(1,length(X));
 
for ii=1:length(X)
  if isstruct(X{ii})
    X{ii}=X{ii}.x;
  end
  [Ti, Ci, ni]=size(X{ii});
  
  if ii==1
    n = ni;
  elseif ni~=n
    error('All blocks should have the same length');
  end
  
  if Ti~=1
    dd1(ii)=Ti; dd2(ii)=Ci;
  else
    Ci=sqrt(Ci);
    dd1(ii)=Ci; dd2(ii)=Ci;
  end
end

Xout = zeros(dd1*dd2',n);
  
ix0=0;
for ii=1:length(X)
  ix = ix0 + (1:dd1(ii)*dd2(ii));
  ix0= ix(end);
  Xout(ix,:) = opt.weight(ii)*reshape(X{ii},[dd1(ii)*dd2(ii),n]);
end


