function [fval,gg,dval,problem,floss,gloss]=objsqfro(ww, problem, opt)

% function evaluation is counted in lossmn
[fval,gg,dval,problem] = feval(problem.loss, ww, problem, opt);

floss=fval;
gloss=gg;

if problem.hasbias
  ww = ww(1:end-1);
  qq = gg(1:end-1);
  gbias = gg(end);
else
  qq = gg;
  gbias = [];
end

ix0=0;
for ii=1:length(problem.lambda)
  sz  = [problem.ns(ii), problem.nc(ii)];
  ix  = ix0+(1:prod(sz));
  ix0 = ix(end);

  fval = fval + 0.5*problem.lambda(ii)*sum(ww(ix).^2);
  gg(ix) = qq(ix) + problem.lambda(ii)*ww(ix);
  dval = dval - 0.5*sum(qq(ix).^2)/problem.lambda(ii);
end


problem.dpen = norm(qq);