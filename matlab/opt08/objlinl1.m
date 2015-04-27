function [fval,gg,dval,problem,floss,gloss]=objlinl1(ww, problem, opt)

% function evaluation is counted in lossmn
[fval,gg,dval,problem,floss,gloss] = feval(problem.loss, ww, problem, opt);

if problem.hasbias
  ww = ww(1:end-1);
  gbias = gg(end);
  gg = gg(1:end-1);
else
  gbias = [];
end

nm      = abs(ww);
dnm     = abs(gg);

ix1 = find(nm>0);
ix0 = find(nm==0);


gg(ix1) = gg(ix1) + problem.lambda*(ww(ix1)./nm(ix1));
gg(ix0) = gg(ix0).*(1-min(problem.lambda./dnm(ix0), 1));



fval = fval + problem.lambda*sum(nm);
gg   = [gg; gbias];
dval = dval - (sum(nm)*(max(dnm))^2)/(4*problem.lambda);

problem.dpen = max(dnm);
