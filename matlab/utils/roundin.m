function yy=roundin(xx)

xr = rangeof(xx);

xr(1) = ceil(xr(1)-1e-12);
xr(2) = floor(xr(2)+1e-12);

yy = round(xx);
yy = max(yy, xr(1));
yy = min(yy, xr(2));
