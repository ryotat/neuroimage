function out= applyClassifier(fv, model, C, idx)
%out= applyClassifier(fv, classy, C, <idx>)

fv= proc_flaten(fv);
applyFcn= getApplyFuncName(model);


if ~exist('idx','var'), 
  out = feval(applyFcn,C,fv.x);
else
  if isnumeric(fv.x),
    out= feval(applyFcn, C, fv.x(:,idx));
  else
    fv = setTrainset(fv,idx);
    out = feval(applyFcn,C,fv.x);
  end
end






% Dec 2008: copied from IDA toolbox. 
% All rights belong to the authors and Fraunhofer FIRST.IDA.
% http://ida.first.fraunhofer.de/homepages/ida/


