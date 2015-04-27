function cpy= copyStruct(src, varargin)
%cpy= copyStruct(src, f1, ...)
%
% cpy is a copy of struct src, except for fields f1, ... 

cpy= [];
names= fieldnames(src);
for fi= 1:length(names),
  if isempty(strmatch(names{fi}, varargin)),
    cpy= setfield(cpy, names{fi}, getfield(src, names{fi}));
  end
end



% Dec 2008: copied from IDA toolbox. 
% All rights belong to the authors and Fraunhofer FIRST.IDA.
% http://ida.first.fraunhofer.de/homepages/ida/


