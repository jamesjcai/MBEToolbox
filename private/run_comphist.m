

sps={'Afu','an','pm','ca','sc'}

for (k=1:length(sps)),
sp=char(sps{k})


	aln_cds=readfasta(['C:\Documents and Settings\jamescai\Desktop\intragenrep\data\cd_',sp,'.fas'],2,1)
	aln_aa=readfasta(['C:\Documents and Settings\jamescai\Desktop\intragenrep\data\aa_',sp,'.fas'],3,1)

	s=aln_cds.seq;
	s2=aln_aa.seq;

	figure;
	comphist(s,1);
	title(sp);

	%set(gcf,'Title',sp)
	%print -dmeta nt_pm.emf
	print (gcf,'-dmeta', ['nt_',sp,'.emf'])

	figure;
	comphist(s,2);
	title([sp,' - Codon composition']);
	%print -dmeta cd_pm.emf
	print (gcf,'-dmeta', ['cd_',sp,'.emf'])

	figure;
	comphist(s2,3);
	title([sp,' - AA composition']);
	%print -dmeta aa_pm.emf
	print (gcf,'-dmeta', ['aa_',sp,'.emf'])

end