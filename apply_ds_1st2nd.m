function out=apply_ds_1st2nd(C, X)

X=shiftdim(X);
[dd,n]=size(X);

C.W1 = C.Wt*C.W1*C.Ws;
for ii=1:length(C.W2)
  C.W2{ii} = C.Ws2{ii}*C.W2{ii}*C.Ws2{ii}';
end

w = pack1st2nd(C.W1,C.W2)';

out = real(w*X);

if isfield(C,'bias')
  out = out + C.bias;
end
