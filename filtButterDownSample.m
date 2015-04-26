function x = filtButterDownSample(x, order, band, fs, n)

x = filtButter(x, order, band, fs);
x = x(1:n:end, :, :);