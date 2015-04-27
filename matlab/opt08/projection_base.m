% [nmp, fval, dval]=projection_base(nm, C)
function [nmp, fval, dval]=projection_base(nm, C)

% Not checked well yet  
if C<0
  ss = sum(nm);
  C = 1/(ss+C)*ones(1,length(nm));
end

if sum(C.*nm)<=1.0
  nmp = nm;
  fval = 0;
  dval = 0;
  return;
end

  
I = ~isinf(C);

if ~any(I)
  nmp=zeros(size(nm));
  fval = 0;
  dval = 0;
  return;
end

while 1
  lambda = (sum(C(I).*nm(I))-1)/sum(C(I).^2);
  M = zeros(size(nm));
  M(I) = (nm(I)-C(I)*lambda);
  
  if lambda<0
    lambda=0;
    break;
  end

  if all(M>=0)
    break;
  end
  
  I = M>0;
  
end

nmp = zeros(size(nm));
nmp(I) = (nm(I) - lambda*C(I));

I=find(nmp~=0);
fval = 0.5*sum((nm-nmp).^2);
dval = -0.5*sum((nm(I)-lambda*C(I)).^2) +0.5*sum(nm.^2)-lambda;


if abs(fval-dval)>1e-10*max(max(fval,dval),1)
  error('Duality gap larger than 1e-10');
end


% fprintf('Duality gap=%g\n',fval-dval);
