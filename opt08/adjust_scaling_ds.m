function D = adjust_scaling_ds(xx, prob)

if prob.hasbias
  xx = xx(1:end-1);
end

W = reshape(xx, [prob.ns, prob.nc]);

[U,S,V]=svd(W,'econ');
D = U*sqrt(S/sum(diag(S)))*U';
% D = chol(U*S*U'/sum(diag(S)));

% D = sqrtm(W*W');
% D = D/trace(D);

