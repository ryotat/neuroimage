function [W, bias, z, status]=lrds_p300(X, Y, lambda, varargin)
% [W, bias, z]=lrds_p300(X, Y, lambda, bQuiet, precision)
% lambda=exp(linspace(log(0.01), log(100), 20));

opt=propertylist2struct(varargin{:});
opt=set_defaults(opt, 'display', 'final',...
                      'precision', 'default');
[R,C,n] = size(X);

Y = (1+reshape(Y, [6, n/6]))/2;

opt.bQuiet = ~strcmp(opt.display,'iter');

cvx_begin sdp
cvx_gp_precision(1e-5)
cvx_quiet(opt.bQuiet);
cvx_precision(opt.precision);
variable W(R,C);
variable Q1(R,R) symmetric;
variable Q2(C,C) symmetric;
variable z(6,n/6);
minimize -sum(sum(Y.*z))+sum(log(sum(exp(z)))) + lambda*(trace(Q1)+trace(Q2));
subject to
for i=1:n
    trace(W'*X(:,:,i))==z(i);
end
[Q1, -W/2; -W'/2, Q2] >= 0;
cvx_end
status = cvx_status;

bias=0;

if isempty(strfind(status,'Solved'))
  fprintf('CVX_STATUS [%s]. All coefficients set to zero.\n', status);
  W    = zeros(R,C);
  bias = 0;
end
