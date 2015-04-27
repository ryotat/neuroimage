function file_summ = gen_file_summ(label_fmt, subjects, lambda)

sbs = connect_str(subjects, '_');
lms = connect_str(foreach(@num2str, num2cell(lambda([1,end]))),'_');

file_summ=sprintf(strrep(label_fmt, '%g', '%s'),sbs,lms);
