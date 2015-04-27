function X = apply_scaling_l1(D, prob, varargin)

if isnumeric(prob)
  nn   = 1;
  ww   = prob;
  prob = varargin{1};
  X = ww;
else
  [dd, nn]=size(prob.X);
  X = prob.X;
end


if prob.hasbias
  X    = X(1:end-1,:);
  bias = X(end,:);
else
  bias = [];
end

X = diag(D)*X;

if prob.hasbias
  X = [X; bias];
end

