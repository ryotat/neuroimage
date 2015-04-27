function str = num2str_c(num, n)

if ~exist('n','var') || isempty(n)
  n = 2;
end

str = cell(size(num));

for ii=1:prod(size(num))
  number = num(ii);
  expo = floor(log10(number));
  frac = number/10^expo;

  if frac>=10*(1-0.5*10^(-n))
    frac = frac/10;
    expo = expo+1;
  end
  
  
  if frac<1+0.5*10^(1-n)
    front = '';
  else
    front = num2str(frac,n);
    if expo~=0
      front = [front '\times'];
    end
  end

  if expo~=0 || isempty(front)
    back = ['10^{' num2str(expo) '}'];
  else
    back = '';
  end
  

  str{ii} = [front back];
end

if prod(size(str))==1
  str = str{1};
end
