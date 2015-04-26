nRows = 16;
nCols = 32;
n = 1000;
nte = 1000;
sd = 0.0;

X = rand(nRows*nCols,n)*2-1;

w0 = randn(nRows, nCols); 
[U,S,V]=svd(w0); sv=diag(S); nc=ceil(0.3*min(nRows,nCols));
w0 = U(:,1:nc)*spdiag(sv(1:nc))*V(:,1:nc)';

z = reshape(w0,[1,nRows*nCols])*X + sd*randn(1,n);

Y = 2*(z>0)-1;
X = reshape(X, [nRows,nCols,n]);


Xte = rand(nRows*nCols,nte)*2-1;
zte = reshape(w0,[1,nRows*nCols])*Xte + sd*randn(1,nte);

Ytrue = 2*(zte>0)-1;


% Z = reshape(Xte*reshape(W,[nRows*nCols,1]) + sd*randn(ncls*n,1), [ncls, n]);
% [mm, ix] = max(Z);



problem.ns=nRows;
problem.nc=nCols;
problem.X = [reshape(X,[nRows*nCols,n]);ones(1,n)];
problem.Y = Y;
problem.lambda = 10;
problem.hasbias = 1;
problem.linesearch = 'linesearch_backtracking';
problem.obj= 'objsqds';
problem.loss='losslr';
problem.dpen = [];
problem.progress = 'progressds';
problem.projection='truncateds';

