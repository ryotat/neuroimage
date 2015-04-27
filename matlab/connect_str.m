function scon = connect_str(list, con)

scon=cell2mat([foreach(inline(sprintf('[x ''%s'']',con)), list(1:end-1)),list(end)]);
