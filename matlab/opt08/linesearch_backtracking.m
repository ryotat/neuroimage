function [ret, xx, fval, gg, dval, problem, step,floss,gloss]...
    =linesearch_backtracking(xx, fval, gg, dd, dginit, step, problem, opt)

floss=0;
gloss=zeros(size(gg));

if dginit>=0
  fprintf('dg=%g is not a descending direction!\n', dginit);
  [fval,gg,dval, problem, floss, gloss]=feval(problem.obj, xx, ...
                                              problem, opt);
  step = 0;
  ret = -1;
  return;
end


xx0 = xx;
f0  = fval;
gg0 = gg;
cc = 0;

% fprintf('finit=%.20f\n',f0);

while cc<opt.max_linesearch
  ftest = f0  + opt.ftol*step*dginit;
  xx    = xx0 + step*dd;
  [fval, gg, dval, problem, floss, gloss]=feval(problem.obj, xx, problem, opt);
  
  if fval<=ftest
    break;
  end
  % fprintf('[%d] step=%g fval=%.20f > ftest=%.20f\n', cc, step, fval, ftest);
  
  step = step/2;
  cc = cc+1;
end

if cc==opt.max_linesearch
  if opt.display>0
    fprintf('Maximum linesearch=%d reached\n', cc);
  end
  xx   = xx0;
  [fval,gg,dval, problem, floss, gloss]=feval(problem.obj, xx, ...
                                              problem, opt);
  step = 0;
  ret = -2;
  return;
end


ret = 0;