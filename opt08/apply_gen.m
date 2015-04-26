function out=apply_gen(C, X)

X=shiftdim(X);
[d,n]=size(X);

if isfield(C, 'Ww') % symmetric positive definite input
  w = reshape(conj(C.Ww*C.W*C.Ww'), [1, d]);
elseif isfield(C, 'Wt') & isfield(C, 'Ws') % rectangular
  w = reshape(conj(C.Wt*C.W*C.Ws), [1, d]);
else
  w = reshape(conj(C.W), [1, d]);
end

out = real(w*X);

if isfield(C,'bias')
  out = out + C.bias;
end
