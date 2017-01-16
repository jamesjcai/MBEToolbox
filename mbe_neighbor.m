function [tree] = mbe_neighbor(aln,outgrno)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


oldpath=pwd;
if (nargin<2), outgrno=1; end
if (isstruct(aln)), seq=aln.seq; else seq=aln; end
[n,m]=size(seq);
if (outgrno>m), error('NOUTGROUP should smaller than M.'); end

cmd = 'mbe_neighbor.m';
dirstr = fileparts(which(cmd));
cd(dirstr);
if (ispc), sep='\'; else sep='/'; end
infile = [dirstr,sep,'infile'];


[D]=chooseDistance(aln);
i_writematrix(D,aln,infile);

try
    neighbor_mex(outgrno);
    tree=i_readouttree;    
catch
   tree=[]; 
end
%if (neighbor_mex(outgrno)~=0),
%      tree=[];
%else
%      tree=i_readouttree;
%end
cd(oldpath);


function [answer] = i_readouttree
	fid = fopen('outtree', 'r');
	answer = fscanf(fid, '%s');    % It has two rows now.
	fclose(fid);
	x=find(answer==';');
	if ~(isempty(x))
		answer=answer(1:x(1));
	end



function i_writematrix(D,aln,filename)

fid=fopen(filename,'w');
n=size(D,1);
fprintf(fid,'  %d\n',n);
   for i=1:n       
     fprintf(fid,'%10s',i_name10(aln.seqnames{i}));   
   for j=1:n       
     fprintf(fid,'%2.6f ',D(i,j));
   end
     fprintf(fid,'\n');      
   end
fclose(fid);
