% p300_sortbycode - Sorts epochs within subtrials according to the
%                   flashing pattern.
%
% Syntax:
%  [epo, I] = p300_sortbycode(epo)
%
% For example:
%  [8 2 10 7 4 3 11 6 12 9 5 1] [12 4 9 ...]
% ->[1 2 3 4 5 6 7 8 9 10 11 12] [1 2 3...]
%  
%
% Ryota Tomioka 2009

function [epo, I] = p300_sortbycode(epo)

[T,C,n]=size(epo.x);

I = zeros(12,n/12);

for i=1:12
  I(i,:) = find(epo.code==i);
end

I = reshape(I, [1, n]);

epo.x = epo.x(:,:,I);
epo.code = epo.code(:,I);

if isfield(epo, 'y') && ~isempty(epo.y)
  epo.y = epo.y(:,I);
end


