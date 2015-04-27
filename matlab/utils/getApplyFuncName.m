function applyFcn= getApplyFuncName(model)
%applyFcn= getApplyFuncName(model)

if isstruct(model),
  func= getFuncParam(model.classy);
else
  func= getFuncParam(model);
end
applyFcn= ['apply_' func];
if ~exist(applyFcn, 'file'),
  applyFcn= 'apply_separatingHyperplane';
end



% Dec 2008: copied from IDA toolbox. 
% All rights belong to the authors and Fraunhofer FIRST.IDA.
% http://ida.first.fraunhofer.de/homepages/ida/


