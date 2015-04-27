function [fval,gg,dval,problem,out]=lossmnogn(ww, problem, opt)

% Y must be provided as
% Y2=sparse(Y+1, 1:n, ones(1,n), ncls, n);

problem.nfeval = problem.nfeval + 1;  
  
[nsnc, nall]   = size(problem.X);
[ncls, n   ]   = size(problem.Y);


out=zeros(ncls,n);
gg =zeros(nsnc,1);

fval = 0;
dval = 0;

for ii=1:n
  ixtr = (ii-1)*ncls;

  for jj=1:ncls
    ixepo = ixtr + jj;

    for kk=1:problem.ns*problem.nc
      out(jj,ii) = out(jj,ii) + problem.X(kk,ixepo)*ww(kk);
    end
  end
  
  outmax = max(out(:,ii));

  sumexp = 0;
  for jj=1:ncls
    ixepo = ixtr + jj;
    sumexp = sumexp + myexp(out(jj,ii) - outmax);
  end

  fval = fval + (-out(find(problem.Y(:,ii)),ii) + outmax + mylog(sumexp));

  for jj=1:ncls
    ixepo = ixtr+jj;
    logpout = out(jj,ii)-outmax-mylog(sumexp);
    pout    = myexp(logpout);
    alpha   = pout - (jj==find(problem.Y(:,ii)));
    gg = gg + alpha*problem.X(:,ixepo);
    
    dval = dval-pout*logpout;
  end
end

