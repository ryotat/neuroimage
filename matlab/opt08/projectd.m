function [xp, fval, dval]=projectd(x, C)

nm = abs(x);

if sum(C.*nm)<=1.0
  xp = x;
  fval = 0;
  dval = 0;
  return;
end

  
I = ones(1,length(x));

while 1
  lambda = (sum(I.*C.*nm)-1)/sum(I.*C.^2);
  M = I.*(nm - C*lambda);
  
  if all(M>=0)
    break;
  end
  
  I = double(M>0);
end

xp = I.*(x - lambda*sign(x).*C);

fval = 0.5*sum((x-xp).^2);
dval = -0.5*sum((xp~=0).*(nm-lambda*C).^2) +0.5*sum(nm.^2)-lambda;


if abs(fval-dval)>1e-10*fval
  error('Duality gap larger than 1e-10');
else
  fprintf('Duality gap=%g\n',fval-dval);
end
