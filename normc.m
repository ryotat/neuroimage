function nm=normc(C,mode,varargin)

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'v',[],'div', 1, 'dosum', 1);

v = opt.v;

if isnumeric(C)
  C.W=C;
end

if isfield(C,'W')
  switch(mode)
   case {'ds' 'fro'}
    nm=svd(C.W)';
   case 'space'
    nm=sqrt(sum(C.W.^2));
   case 'time'
    nm=sqrt(sum(C.W.^2,2))';
   case 'fourier'
    nm=reshape(abs(C.W),[1,prod(size(C.W))]);
   otherwise
    nm=norm(C.W,mode);
  end
elseif isfield(C,'W1')
  nm=[];
  W=[{C.W1}, C.W2];
  if isempty(v)
    v=ones(1,length(W));
  end
  if opt.div==1
    v=v/v(1);
  end
  switch(mode)
   case {'ds' 'fro'}
    for ii=1:length(W);
      nm = [nm,sqrt(v(ii))*svd(W{ii})'];
    end
   otherwise 
    error('mode: %s not supported.',mode);
% $$$     
% $$$    case 'fro'
% $$$     for ii=1:length(W);
% $$$       nm = nm+v(ii)*norm(W{ii},'fro')^2;
% $$$     end
% $$$     nm=sqrt(nm);
  end
end

if opt. dosum
  switch(mode)
   case {'ds', 'space', 'time'}
    nm = sum(nm);
   case 'fro'
    nm = sqrt(sum(nm.^2));
   otherwise 
    error('mode: %s not supported.',mode);
  end
end
