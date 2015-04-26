function [fval,gg,dval,problem]=losslr(ww, problem, opt)

problem.nfeval = problem.nfeval + 1;  
  
[nsnc, n]   = size(problem.X);


z     = (ww'*problem.X).*problem.Y;
pout  = zeros(size(z));

I1 = find(z>=0);
I2 = find(z<0);

pout(I1)= 1./(1+exp(-z(I1)));
pout(I2)= exp(z(I2))./(1+exp(z(I2)));

IP = find(0 < pout & pout < 1);

alpha   = (-problem.Y.*(1-pout))';

fval    = sum(log(1+exp(-z(I1)))) + sum(-z(I2)+log(exp(z(I2))+1));
gg      = problem.X*alpha;
dval    = -sum(pout(IP).*log(pout(IP)) + (1-pout(IP)).*log(1-pout(IP)));
