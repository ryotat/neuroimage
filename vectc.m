function xx = vectc(C)

if isfield(C,'W')
  xx = reshape(C.W,[prod(size(C.W)),1]);
elseif isfield(C,'W1')
  xx = reshape(C.W1,[prod(size(C.W1)),1]);
  for ii=1:length(C.W2)
    xx = [xx; reshape(C.W2{ii},[prod(size(C.W2{ii})),1])];
  end
end

if isfield(C,'bias')
  xx = [xx; C.bias];
end
