function [fval,gg,dval,problem,floss,gloss]=lossmn(ww, problem, opt)

% Y must be provided as
% Y2=sparse(Y+1, 1:n, ones(1,n), ncls, n);

problem.nfeval = problem.nfeval + 1;  
  
[nsnc, nall]   = size(problem.X);
[ncls, n   ]   = size(problem.Y);

out     = reshape(ww'*problem.X, [ncls, n]);
outmax  = max(out);
sumexp  = sum(exp(out-ones(ncls,1)*outmax));

logpout = out-ones(ncls,1)*(outmax+log(sumexp));
pout    = exp(logpout);

alpha   = reshape(-problem.Y + pout, [nall,1]);


fval    = sum(-sum(problem.Y.*out)+outmax+log(sumexp));
gg      = problem.X*alpha;
dval    = -sum(sum(pout.*logpout));

floss = fval;
gloss = gg;
