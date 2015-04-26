function ww = truncatel1(ww, problem, opt)
  
if problem.hasbias
  bias = ww(end);
  ww=ww(1:end-1);
else
  bias = [];
end

nm  = abs(ww);
ix0 = find(nm<=opt.epsh);

if length(ix0)>0
  fprintf('Zeroing %d components:', length(ix0));
end


ww(ix0)=0;


ww = [ww; bias];
