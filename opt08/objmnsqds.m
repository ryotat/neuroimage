function [fval,gg,dval,problem]=objmnsqds(ww, problem, opt)

% function evaluation is counted in lossmn
[fval,gg,dval,problem] = lossmn(ww, problem, opt);

if problem.hasbias
  ww = ww(1:end-1);
  Q = reshape(gg(1:end-1), [problem.ns, problem.nc]);
  gbias = gg(end);
else
  Q = reshape(gg, [problem.ns, problem.nc]);
  gbias = [];
end
ww = reshape(ww, [problem.ns, problem.nc]);

[U,S,V] = svd(ww);
nm      = diag(S);
dnm     = svd(Q);

lmnm = problem.lambda*sum(nm);

ix1 = find(nm>0);
ix0 = find(nm==0);

U1=U(:,ix1); V1=V(:,ix1);
U0=U(:,ix0); V0=V(:,ix0);


if isempty(ix0)
  Q = Q + lmnm*U1*V1';
else
  [Uq,Sq,Vq] = svd(U0'*Q*V0);
  sq = diag(Sq);

  if lmnm>0
    cc = max(-1, -sq/lmnm);
  else
    cc = -ones(1,length(sq));
  end
  Q = Q + lmnm*(U1*V1' + U0*Uq*spdiag(cc)*Vq'*V0');
end




fval = fval + 0.5*problem.lambda*(sum(nm))^2;
gg   = [reshape(Q, [problem.ns*problem.nc,1]); gbias];
dval = dval - 0.5*(max(dnm))^2/problem.lambda;
