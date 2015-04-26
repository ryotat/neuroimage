function [xx, gg]=readxg(file)
  
fid=fopen(file);

xx = sscanf(fgets(fid), '%f ');
gg = sscanf(fgets(fid), '%f ');

fclose(fid);