function N = countComp(S, varargin)

opt = propertylist2struct(varargin{:});
opt = set_defaults(opt, 'tolabs', 0, ...
                        'tolrel', 0.01,...
                        'dosum', 1);

S = abs(S);

if size(S,1)==1
  S=S';
end

N=zeros(size(S,1),size(S,2));

M=max(S,[],1);

I=M>opt.tolabs;
N(:,I)=S(:,I)>opt.tolrel*ones(size(S,1),1)*M(I);
N(:,~I)=0;

if opt.dosum
  N = sum(N,1);
end


