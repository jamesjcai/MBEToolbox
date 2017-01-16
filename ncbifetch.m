function ncbifetch(accession,db)
%NCBIFETCH - fetches CDS sequences from NCBI database for given genbank accessions
%
%[s,aln] = ncbifetch(accession,db)
%
% e.g., ncbifetch({'47824786','71001358'},'pep2cds')
%%

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



if (nargin<2)
      db='pep2cds';
      defaccession={'70997565,67515801,39975839'};
	prompt = {'Enter GenPept Accession(s):'};
	dlg_title = 'Fetch CDS sequences from NCBI';
	num_lines = 1;
	answer = inputdlg(prompt,dlg_title,num_lines,defaccession);
	if (isempty(answer)),    return; end

	x=strread(char(answer{1}),'%d','delimiter',',');
	%x=unique(x);
	accession=num2cell(x);
end

if ~(iscell(accession)), accession={accession}; end

base='http://www.ncbi.nlm.nih.gov';
n=length(accession);

s='';
for (k=1:n),
	acc=num2str(accession{k});
	switch (lower(db))
	    case ('pubmed')
		 url = strcat(base, '/entrez/query.fcgi?cmd=Text&db=PubMed&dopt=MEDLINE&uid=', acc);
		 [s0,status]=urlread(url);
		 s0=i_removepre(s);

	    case ('pep2cds')
		url = strcat(base, '/entrez/viewer.fcgi?val=',acc);
		disp(sprintf('... getting CDS accession from %s',url))

		[s0,status]=urlread(url);
		[cdsacc,cdsitem] = i_getcdsacc2(s0);
		if ~(isempty(cdsacc))
			%newurl = strcat(base,'/entrez/viewer.fcgi?db=nucleotide&qty=1&c_start=1&list_uids=',cdsacc,'&itemID=1&dopt=fasta&dispmax=5&sendto=t');
			newurl = strcat(base,'/entrez/viewer.fcgi?db=nucleotide&qty=1&c_start=1&list_uids=',cdsacc,'&itemID=',cdsitem,'&dopt=fasta&dispmax=5&sendto=t');
			disp(sprintf('...... fetching from %s',newurl))
			[s0,status]=urlread(newurl);
		else
			disp(sprintf('...... cannot parse CDS information for %s',acc))
			%return
			s0=' ';
		end
	    otherwise
		url = strcat(base, '/entrez/viewer.fcgi?val=',acc);
		[s0,status]=urlread(url);
		[cdsacc] = i_getcdsacc(s0);
        newurl = strcat(base,'/entrez/viewer.fcgi?db=nucleotide&qty=1&c_start=1&list_uids=',cdsacc,'&itemID=1&dopt=fasta&dispmax=5&sendto=t');
		[s0,status]=urlread(newurl);
	end

	if (s0(end)>65),   %not return LF(10)/CR(13) then add an return
		s0=sprintf('%s\n.',s0);
	end
	s=[s,s0];
end


    [filename, pathname,filterindex] = uiputfile( ...
       {'*.fasta;*.fas', 'FASTA Format Files (*.fasta, *.fas)';
        '*.*',  'All Files (*.*)'}, ...
        'Save as');
	if ~(filename), return; end
	filename=[pathname,filename];
	if (filterindex==1)
		if (isempty(find(filename=='.')))
		filename=[filename,'.fas'];
		end
	end

       [fid,Msg] = fopen(filename,'wt');
       if fid == -1, error(Msg); end
       fprintf(fid,'%s',s);
       fclose(fid);





% 	if (nargout>1)   % need Aln
% 	       tempf = tempfilename;
% 	       [fid,Msg] = fopen(tempf,'wt');
% 	       if fid == -1, error(Msg); end
% 	       fprintf(fid,'%s',s);
% 	       fclose(fid);
%
% 	       if (n==1),
% 	          aln=readfasta(tempf,2,1);
% 	       elseif (n>1)
% 	       	   ButtonName=questdlg('Do you want to align downloaded CDSs?', ...
% 				    'CDS sequences downloaed', ...
% 				    'Yes','Cancel','Yes');
% 		    switch ButtonName,
% 		    case 'Yes',
% 			aln = alignseqfile('CDS',tempf,1);
% 		    otherwise
% 		        aln=[];
% 		    end
% 		end
% 	end





function [acc] = i_getcdsacc(s)
acc='';
x=strread(s,'%s','delimiter','\n');
for (k=1:length(x)),
      linex=char(x(k));
      %fprintf('%s\n',linex);
      if ~(isempty(findstr(linex,'CDS'))),
	      % '     <a href=/entrez/viewer.fcgi?val=67539623&itemID=1&view=gbwithparts>CDS</a>             1..2371'
	      % fprintf('%s\n',linex);
	      [start,finish,tokens] = regexp(linex,'[0-9]+');
	      acc=linex([start(1):finish(1)]);

	      return
      end

end


function [acc,item] = i_getcdsacc2(s)
acc='';
item='';
x=strread(s,'%s','delimiter','\n');
for (k=1:length(x)),
      linex=char(x(k));
      %fprintf('%s\n',linex);
      if ~(isempty(findstr(linex,'CDS'))),
	      % '     <a href=/entrez/viewer.fcgi?val=67539623&itemID=1&view=gbwithparts>CDS</a>             1..2371'
	      % fprintf('%s\n',linex);
	      [start,finish,tokens] = regexp(linex,'[0-9]+');
	      acc=linex([start(1):finish(1)]);
	      item=linex([start(2):finish(2)]);
	      return
      end

end



function [s] = i_removepre(s)
s = strrep(s,'<pre>','');
s = strrep(s,'</pre>','');
