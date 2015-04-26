function r = rangeof(X)

n = prod(size(X));
X = reshape(X, [n, 1]);

r = [min(X), max(X)];