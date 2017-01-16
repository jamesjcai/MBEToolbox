function showseq61(s61,flag)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

if (nargin<2)
      flag=0;
end

	codon={'AAA' 'AAC' 'AAG' 'AAT' 'ACA' 'ACC' 'ACG' 'ACT' 'AGA' 'AGC' 'AGG' 'AGT' 'ATA' 'ATC' 'ATG' 'ATT'...
	      'CAA' 'CAC' 'CAG' 'CAT' 'CCA' 'CCC' 'CCG' 'CCT' 'CGA' 'CGC' 'CGG' 'CGT' 'CTA' 'CTC' 'CTG' 'CTT'...
	      'GAA' 'GAC' 'GAG' 'GAT' 'GCA' 'GCC' 'GCG' 'GCT' 'GGA' 'GGC' 'GGG' 'GGT' 'GTA' 'GTC' 'GTG' 'GTT'...
	      'TAA' 'TAC' 'TAG' 'TAT' 'TCA' 'TCC' 'TCG' 'TCT' 'TGA' 'TGC' 'TGG' 'TGT' 'TTA' 'TTC' 'TTG' 'TTT'};


icode=1;
[TABLE,CODON] = codontable;
stops=find(TABLE(icode,:)=='*');
codon(stops)=[];
TAB=TABLE(icode,:);

[n,m]=size(s61);


switch (flag)
    case (0)
	for (i=1:n),
		fprintf('seq%4d: ',i);
	for (j=1:m),
		id=s61(i,j);
		c=char(codon(id));
		a=translate(c);
		fprintf('%s (%s) ', c,a);
	end
		fprintf('\n');
	end
    case (1)
	for (i=1:m),
		fprintf('site%4d: ',i);
	for (j=1:n),
		id=s61(j,i);
		c=char(codon(id));
		a=translate(c);
		fprintf('%s (%s) ', c,a);
	end
		fprintf('\n');
	end
   otherwise 
	 error('flag error')
end
