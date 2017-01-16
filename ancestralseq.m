function [seqout]=ancestralseq(aln,usecodeml)
%ANCESTRALSEQ - reconstructs ancestral sequence with PAML.
%
%[seqout]=ancestralseq(aln,method)
% method=1, BASEML; method=2, CODEML.
%
%Sequence type of SEQOUT will be the same as ALN.SEQTYPE.

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


if (nargin<2), usecodeml=1; end
if size(aln.seq,1)>2
    error('This function is only working for alignment with two sequences.');
end

seqout=[];
oldpath=pwd;
cdmbe;
cd 'addins';
cd 'paml';

% dir
if size(aln.seq,1)>2
    [tree] = gettreedlg(aln);
else
    tree='(1, 2);';
end
fid=fopen('treefile','w');
fprintf(fid,'%s\n',tree);
fclose(fid);
writephylip_s(aln,'infile');

if (isunix), slas='./'; else slas=''; end
if (usecodeml==1)
	cmd=[slas,'codeml'];
	if (aln.seqtype)~=2, error('Needs codon sequence.'); end
else
	cmd=[slas,'baseml'];
end

try     
    disp(['Running command: ',cmd])
	[s,w] = system(cmd);
    disp(w);
catch
    return;
end


txt = textread('rst','%s','delimiter','\n','whitespace','');
idx=[];
for k=1:length(txt)
if strfind(txt{k},'node #')==1
    %seq=txt{k};
    idx=[idx;k];
end
end
if ~isempty(idx)
    seq=txt{idx(1)};
else
    seq='';
end

idxnum=regexp(seq,'\d');
%seq2=strtrim(seq(max(idxnum)+1:end));
seq2=seq(max(idxnum)+1:end);
seq2=seq2(seq2~=' ');

if aln.seqtype==3
    seqout=i_encode_a(seq2);    
else
    seqout=i_encode_n(seq2);
end
cd(oldpath);
