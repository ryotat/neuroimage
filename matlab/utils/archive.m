function S=archive(varargin)

S = [];
for i=1:length(varargin)
  name =varargin{i};
  S = setfield(S, name, evalin('caller', name));
end
