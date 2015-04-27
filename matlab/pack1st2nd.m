function X=pack1st2nd(X1, X2)

  fprintf('In pack1st2nd length(X2)=%d\n', length(X2));  

if isstruct(X1)
  X1=X1.x;
end
 
[T,C,n]=size(X1);

X = zeros(T*C+length(X2)*C^2,n);
X(1:T*C,:)=reshape(X1, [T*C,n]);
  
ix0=T*C;
for ii=1:length(X2)
  ix = ix0 + (1:C^2);
  ix0= ix(end);
  
  if isstruct(X2{ii})
    X(ix,:) = reshape(X2{ii}.x,[C*C,n]);
  else
    X(ix,:) = reshape(X2{ii},[C*C,n]);
  end
end


