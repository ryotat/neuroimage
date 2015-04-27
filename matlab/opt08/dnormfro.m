function dval = dnormfro(ww, gg, problem, opt)
if problem.hasbias
  gg = gg(1:end-1);
end

len = problem.ns.*problem.nc;

ix1= [1, 1+cumsum(len(1:end-1))];
ix2= cumsum(len);

nm2 = zeros(1,length(len));
for ii=1:length(len)
  nm2(ii) = sum(gg(ix1(ii):ix2(ii)).^2);
end

dval = -sqrt(sum(nm2./(problem.C.^2)));
  