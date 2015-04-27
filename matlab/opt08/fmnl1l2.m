function [f, g, fl, gl] = fmnl1l2(W, X, Y, lambda, step)

if ~exist('step','var') 
  step = 0;
end

  
[ns, nc, nall] = size(X);
n              = size(Y,2);
ncls           = nall/n;

fl = 0;
gl = zeros(ns,nc);
fn = 0;
gn = zeros(ns,nc);

Y2=zeros(ncls,n);
for i=1:n
  Y2(Y(i)+1,i)=1;
end


out = reshape(reshape(W, [1,ns*nc])*reshape(X, [ns*nc, nall]), [ncls, n]);
outmax = max(out);

sumexp = sum(exp(out-ones(ncls,1)*outmax));

fl = sum(-sum(Y2.*out)+outmax+log(sumexp));

coeff = reshape(-Y2 + exp(out-ones(ncls,1)*(outmax+log(sumexp))), [nall,1]);

gl = reshape(reshape(X,[ns*nc,nall])*coeff,[ns,nc]);

for ii=1:nc
  nm  = norm(W(:,ii));
  nmg = norm(gl(:,ii));
  
  fn = fn + lambda*nm;
  
  if step==0 && nm<1e-6 || step>0 && nm==0
    if nmg <= lambda
      gn(:,ii)=-gl(:,ii);
    else
      gn(:,ii)=-lambda/nmg*gl(:,ii);
    end
  else
    gn(:,ii) = lambda/nm*W(:,ii);
  end
end


f = fl+fn;
g = gl+gn;
