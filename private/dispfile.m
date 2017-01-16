function [str]=dispfile(textfile)

fid=fopen(textfile,'r');
if (fid==-1) 
  error('Please check your filename, cannot open file');
end

id=1;
while 1
     tline = fgetl(fid);
     if ~ischar(tline), break, end
     disp(tline)
     if (nargout>0), str{id}=tline; id=id+1; end
end
fclose(fid);