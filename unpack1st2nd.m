function [X1,X2]=unpack1st2nd(X,T,C,n,varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'flat', 0,'weight',[]);

len = (size(X,1)-T*C)/(C^2);

if isempty(opt.weight)
  opt.weight = ones(1,len+1);
elseif length(opt.weight)~=len+1
  error('Opt.weight must have the same length as X.');
end

X1 = opt.weight(1)*reshape(X(1:T*C,:),[T,C,n]);


X2 = cell(1,len);

ix0=T*C;
for ii=1:len
  ix =ix0+(1:C^2);
  ix0=ix(end);
  if opt.flat
    X2{ii} = opt.weight(ii+1)*reshape(X(ix,:),[1,C*C,n]);
  else
    X2{ii} = opt.weight(ii+1)*reshape(X(ix,:),[C,C,n]);
  end
end



