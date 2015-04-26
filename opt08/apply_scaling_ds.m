function X = apply_scaling_ds(D, prob, varargin)

if isnumeric(prob)
  nn   = 1;
  ww   = prob;
  prob = varargin{1};
  
  if prob.hasbias
    X    = ww(1:end-1);
    bias = ww(end);
  else
    X    = ww;
    bias = [];
  end
  
else
  [dd, nn]=size(prob.X);
  if prob.hasbias
    X    = prob.X(1:end-1,:);
    bias = 1;
  else
    X    = prob.X;
    bias = [];
  end
  
end

X = reshape(X, [prob.ns, prob.nc*nn]);
X = reshape(D*X, [prob.ns*prob.nc,nn]);

if prob.hasbias
  X = [X; bias*ones(1,nn)];
end

