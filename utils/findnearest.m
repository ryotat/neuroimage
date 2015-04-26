function ix=findnearest(x, xi)
% findnearest - find elements in x nearest to xi
%
% ix=findnearest(x, xi)

if isnan(xi)
  ix = [];
else
  ix = interp1(x, 1:length(x), xi, 'nearest','extrap');
end
