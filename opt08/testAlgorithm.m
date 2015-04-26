function [result, nm] = testAlgorithm(nList, lambdaList, problem)

if isfield(problem, 'lambda')
  fldname = 'lambda';
else
  fldname = 'C';
end

if isempty(problem.ncls)
  ncls = 1;
else
  ncls = problem.ncls;
end

if ~isempty(strfind(problem.suf, 'con'))
  lambda_exp = -1;
else
  lambda_exp = 1;
end


X = problem.X;
Y = problem.Y;


x0 = zeros(size(X,1),1);

for ii=1:length(nList)
  n = nList(ii);
  ix = 1:ncls*n;
  problem.X = X(:,ix);
  problem.Y = Y(:,1:n);
  
  for jj=1:length(lambdaList)
    lambda=lambdaList(jj);
    problem=setfield(problem, fldname, lambda^lambda_exp);
    file = sprintf('logfile_%d_%g_%s.txt', n, lambda, problem.suf);
    
    fprintf('[n=%d lambda=%g file=%s]\n',n,lambda,file);
    
    [xx,status]=feval(problem.solver,x0,problem,'display',2,'logfile', ...
                      file,'maxiter',10000);
    W=reshape(xx,[problem.ns,problem.nc]);
    if ~isempty(strfind(problem.suf,'gl'))
      nm=sqrt(sum(W.^2));
    elseif ~isempty(strfind(problem.suf,'ds'))
      nm=svd(W);
    else
      nm=[];
    end
    result(ii,jj)=archive('n','lambda','W','status','nm');
    
    fprintf('\nsum(nm)=%g\n\n',sum(nm));
  end
end

nm=getfieldarray(result,'nm');

  
