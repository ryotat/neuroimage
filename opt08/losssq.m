function [fval,gg,dval,problem,floss,gloss]=losssq(ww, problem, opt)

problem.nfeval = problem.nfeval + 1;  
  
[dd, n]   = size(problem.X);


alpha  = (ww'*problem.X-problem.Y)';

fval   = sum(alpha.^2)/2;
gg     = problem.X*alpha;
dval   = -sum(alpha.^2)/2-problem.Y*alpha;

if problem.hasbias
  dval = dval + ww(end)*sum(alpha);
  % fprintf('* bias*sum(alpha)=%g\n',ww(end)*sum(alpha));
end

floss  = fval;
gloss  = gg; 