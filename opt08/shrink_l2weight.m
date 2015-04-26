function [yy, fval, dval]=shrink_l2weight(xx, C, len, print)

if ~exist('len','var') || isempty(len)
  len = length(xx);
end

if ~exist('print','var') || isempty(print)
  print=0;
end

if length(C)~=length(len)
  error('The length of C and len must be the same.');
end


C = C.^2;

ix1= [1, 1+cumsum(len(1:end-1))];
ix2= cumsum(len);
nm2 = zeros(1,length(len));
for ii=1:length(len)
  nm2(ii) = sum(xx(ix1(ii):ix2(ii)).^2);
end

if sum(nm2.*C)<=1.0
  yy = xx;
  fval=0;
  dval=0;
  return;
end

% keyboard;
lm = norm(xx)^2/4;
yy=zeros(size(xx));
count = 1;
while 1
  nmy = 0;
  for ii=1:length(len)
    yy(ix1(ii):ix2(ii)) = xx(ix1(ii):ix2(ii))/(1+lm*C(ii));
    nmy = nmy + C(ii)*sum(yy(ix1(ii):ix2(ii)).^2);
  end
  fval = 0.5*sum((xx-yy).^2);
  dval = 0.5*sum(lm*C.*nm2./(1+lm*C))-0.5*lm;
  g = 0.5*sum(C.*nm2./((1+lm*C).^2))-0.5;
  h = -sum((C.^2).*nm2./((1+lm*C).^3));

  if print>1
    fprintf('***[%d] lm=%g fval=%g dval=%g |g|=%g nmy=%g\n',count,lm,fval,dval, ...
            norm(g), nmy);
  end
  
  if norm(g)<1e-6 && nmy<=1.0 || count>100
    break;
  end
  
  step=1.0;
  while 1
    lm1 = lm -step*g/h;
    dval1 = 0.5*sum(lm1*C.*nm2./(1+lm1*C))-0.5*lm1;
    if dval1>=dval && lm1>=0
      lm = lm1;
      break;
    end
    step = step/2;
  end
    
  count = count+1;
end

if print>0
  fprintf('***[%d] lm=%g fval=%g dval=%g |g|=%g nmy=%g\n',count,lm,fval,dval, norm(g), nmy);
end
