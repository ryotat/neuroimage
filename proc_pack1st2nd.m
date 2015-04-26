function fv=proc_pack1st2nd(fv, varargin)

fprintf('In proc_pack1st2nd length(varargin)=%d\n', length(varargin));  

[T,C,n]=size(fv.x);

if length(varargin)>0
  X = pack1st2nd(fv.x, varargin);
else
  fvcov = proc_covariance(fv);
  X = [reshape(fv.x, [T*C,n]); shiftdim(fvcov.x)];
end


fv.x = X;

fv.Nsamples = T;
fv.Nchannels = C;
