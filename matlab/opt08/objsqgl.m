function [fval,gg,dval,problem,floss,gloss]=objsqgl(ww, problem, opt)

% function evaluation is counted in lossmn
[fval,gg,dval,problem,floss,gloss] = feval(problem.loss, ww, problem, opt);

if problem.hasbias
  ww = ww(1:end-1);
  gbias = gg(end);
  gg = gg(1:end-1);
else
  gbias = [];
end

ww = reshape(ww, [problem.ns, problem.nc]);
gg = reshape(gg, [problem.ns, problem.nc]);

nm      = sqrt(sum(ww.^2));
dnm     = sqrt(sum(gg.^2));

lmnm = problem.lambda*sum(nm);

ix1 = find(nm>0);
ix0 = find(nm==0);


gg(:,ix1) = gg(:,ix1) + lmnm*(ww(:,ix1)/spdiag(nm(ix1)));
gg(:,ix0) = gg(:,ix0)*spdiag(1-min(lmnm./dnm(ix0), 1));



fval = fval + 0.5*problem.lambda*(sum(nm))^2;
gg   = [reshape(gg, [problem.ns*problem.nc,1]); gbias];
dval = dval - 0.5*(max(dnm))^2/problem.lambda;

problem.dpen = max(dnm);