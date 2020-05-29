function n = linecount(fid)
% linecount(fid)
% count lines in the file with fileID fid
n = 0;
tline = fgetl(fid);
while ischar(tline)
  tline = fgetl(fid);
  n = n+1;
end