function ww = projl1(ww, problem, opt)
nm = abs(ww)';

nmp = projection_base(nm, problem.C)';

ix = nmp>0;

ww(ix)  = ww(ix).*nmp(ix)./nm(ix)';
ww(~ix) = 0;
