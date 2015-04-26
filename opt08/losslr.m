function [fval,gg,dval,problem,floss,gloss]=losslr(ww, problem, opt)

problem.nfeval = problem.nfeval + 1;  
  
[nsnc, n]   = size(problem.X);


z       = (ww'*problem.X).*problem.Y;
z2      = 0.5*[z; -z];
outmax  = max(z2);
sumexp  = sum(exp(z2-ones(2,1)*outmax));
logpout = z2-ones(2,1)*(outmax+log(sumexp));;
pout    = exp(logpout);


alpha   = (-problem.Y.*pout(2,:))';

fval    = -sum(logpout(1,:));
gg      = problem.X*alpha;
dval    = -sum(sum(pout.*logpout));

if problem.hasbias
  dval = dval + ww(end)*sum(alpha);
  % fprintf('* bias*sum(alpha)=%g\n',ww(end)*sum(alpha));
end

floss = fval;
gloss = gg;

