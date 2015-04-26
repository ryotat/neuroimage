function [fval,gg,dval,problem,floss,gloss]=objlinds(ww, problem, opt)

% function evaluation is counted in lossmn
[fval,gg,dval,problem] = feval(problem.loss, ww, problem, opt);

floss=fval;
gloss=gg;

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

ix1 = find(nm>opt.epssvd);
ix0 = find(nm<=opt.epssvd);

U1=U(:,ix1); V1=V(:,ix1);
U0=U(:,ix0); V0=V(:,ix0);

if isempty(ix0)
  Q = Q + problem.lambda*U1*V1';
else
  if problem.ns<=problem.nc
    ix1c=setdiff(1:problem.nc, ix1);
    V1c=V(:,ix1c);
    [Uq,sq,Vq] = svd(U0'*Q*V1c);
    if length(ix0)>1
      sq = diag(sq);
    end
    
    cc = -min(problem.lambda, sq);

    Q = Q + (problem.lambda*U1*V1' + U0*Uq*spdiags(cc,0,length(ix0),length(ix1c))*Vq'*V1c');
  else
    ix1c=setdiff(1:problem.ns, ix1);
    U1c=U(:,ix1c);
    [Uq,sq,Vq] = svd(U1c'*Q*V0);
    if length(ix0)>1
      sq = diag(sq);
    end
    
    cc = -min(problem.lambda, sq);

    Q = Q + (problem.lambda*U1*V1' + U1c*Uq*spdiags(cc,0,length(ix1c),length(ix0))*Vq'*V0');
  end
  % subplot(1,2,1); imagesc(U'*Q*V); colorbar;
  % subplot(1,2,2); plot(diag(U'*Q*V));
end



fval = fval + problem.lambda*sum(nm);
gg   = [reshape(Q, [problem.ns*problem.nc,1]); gbias];
dval = dval - sum(nm)*(max(dnm))^2/(4*problem.lambda);

problem.dpen = max(dnm);