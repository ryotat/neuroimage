% mnl1l2 - Solves group penalized (L1-L2) multinomial regression.
%          The feature vector should consist of nc blocks of equal
%          size (ns). The likelihood is defined as the softmax function:
% 
%          P(y=k|x) = exp(<w,x_k>)/sum_{k'}(exp(<w,x_k'>)) (k=0,...,ncls-1)
%
%          The problem is to solve the following optimization:
%
%          minimize sum_i(-log(P(y=k^i|x^i))) + lambda * sum(sqrt(sum(w.^2)))
%       
% Syntax:
%   [w, z, status] = mnl1l2(X, Y, lambda, <w0>, <lsmethod>, <display>)
%
% Inputs:
%   X      :  input vectors.                [ns, nc, ncls*ntrials]
%   Y      :  class labels (0,...,ncls-1).  [1, ntrials]
% lambda   :  regularization constant
%   w0     :  weight vectors (optional).    [ns, nc]
% lsmethod :  line-search method (optional; default=1)
%             0 = LBFGS_LINESEARCH_MORETHUENTE
%             1 = LBFGS_LINESEARCH_BACKTRACKING
% display  :  0 = no display
%             1 = final
%             2 = iteration
%
% Outputs:
%   w      : weight vectors.                [ns, nc]
%   z      : classifier outputs.            [1, ncls*ntrials]
% status   : miscellaneous information.     [struct]
%
%
% To compile:
% >> mex mnl1l2.c lbfgs.c
%  (requires lbfgs.h and arithmetic_ansi.h)
%
% See also:
% mnl1l2.c, lbfgs.c train_mnl1l2.m, apply_mnl1l2.m
%
% Ryota Tomioka 2008