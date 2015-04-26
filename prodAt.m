function Y = prodAt(X, A, dim)
% prodAt - multiply A at the dim'th dimension of a tensor X
% 
% Y = prodAt(X, A, dim)
%
nd = ndims(X);

if dim~=1
  X = permute(X, [dim, 1:dim-1, dim+1:nd]);
end
sz = size(X);
Y = reshape(A*reshape(X,[sz(1),prod(sz(2:end))]), [size(A,1), sz(2:end)]);

if dim~=1
  Y = permute(Y, [2:dim, 1, dim+1:nd]);
end
