% p300_decode - Decode Farwell & Donchin speller system.
%               assumes out is organized into blocks of 12 flashing
%               patterns, namely:
%               [1st col,...,  6th col, 1st row, ... 6th row] [1st col,...]...
%               Each output value is the log likelihood of selecting
%               the corresponding row or column.
% Syntax:
%  decode = p300_decode(out,Tab,sign)  
%
% Ryota Tomioka 2009

function decode = p300_decode(out,Tab,sign)
  
if ~exist('sign','var')
  sign = 1;
end

sz = size(shiftdim(out));

out = reshape(out, [12, prod(sz)/12]);

[mm, ixc] = max(sign*out(1:6,:));
[mm, ixr] = max(sign*out(7:12,:));

decode= diag(Tab(ixr,ixc))';

if length(sz)>2 || sz(1)>12
  sznew = sz;
  if sznew(1)>12
    sznew(1) = sznew(1)/12;
  else
    sznew = sznew(2:end);
  end
  decode= reshape(decode, sznew);
end


