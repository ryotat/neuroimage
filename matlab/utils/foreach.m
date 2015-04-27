function B=foreach(fun, A, depth, varargin)
% foreach - apply function to all elements in a cell array
% B=foreach(fun, A, depth, varargin)
%  depth : max depth of evaluation for nested cell matrices
%          (default inf)


if ~exist('depth','var') | isempty(depth)
  depth = inf;
end

n=prod(size(A));
B=cell(size(A));

for i=1:n
  if ~iscell(A{i}) | depth == 0
    B{i} = feval(fun, A{i}, varargin{:});
  else
    B{i} = foreach(fun, A{i}, depth-1, varargin{:});
  end
end