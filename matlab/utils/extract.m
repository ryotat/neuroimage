function extract(S)

names = fieldnames(S);

for i=1:length(names)
  assignin('caller', names{i}, getfield(S, names{i}));
end
