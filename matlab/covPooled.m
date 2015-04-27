function C = covPooled(X)
%% C = covPooled(X)
%%
%% Input:
%%   X : T*m*n with T samples m variables and n conditions


[T,m,n] = size(X);


C = zeros(m);
for i=1:n
  C = C + cov(X(:,:,i));
end

C = C/n;