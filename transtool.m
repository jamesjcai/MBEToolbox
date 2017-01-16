function transtool()

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

   prompt={'FASTA format:'};
   def={''};
   dlgTitle='Sequence Input';
   lineNo=10;
   AddOpts.Resize='off';
   AddOpts.WindowStyle='modal';
   AddOpts.Interpreter='none';
   answer=inputdlg(prompt,dlgTitle,lineNo,def,AddOpts);
   
 if ~(isempty(answer)),     
      seq=cellstr(answer{1});
      n=length(seq);
     if (n>=2&length(seq(2))>0)
       tempf = tempfilename;
       [fid,Msg] = fopen(tempf,'wt');
       if fid == -1, error(Msg); end
       for (k=1:n),
           fprintf(fid,'%s\n',char(seq(k)));
       end
       fclose(fid);

       aln=readfasta(tempf,2,1);
       aln2=translatealn(aln);
       viewseq(aln2)
end
end