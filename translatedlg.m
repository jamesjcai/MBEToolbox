function [alndna,alnpro] = translatedlg

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

   prompt={'FASTA format:'};
   def={''};
   dlgTitle='Sequence Translation Dialog Box';
   lineNo=10;
   AddOpts.Resize='off';
   AddOpts.WindowStyle='modal';
   AddOpts.Interpreter='none';
   answer=inputdlg(prompt,dlgTitle,lineNo,def,AddOpts);
   
 if ~(isempty(answer)),     
      seq=cellstr(answer{1});
      n=length(seq);
     if (n>=2&&~isempty(seq(2)))
       if (answer{1}(1)~='>'), error('Sequence is not in FASTA format'); end           
       tempf = tempfilename;
       [fid,Msg] = fopen(tempf,'wt');
       if fid == -1, error(Msg); end
       for (k=1:n),
           fprintf(fid,'%s\n',char(seq(k)));
       end
       fclose(fid);

       alndna=readfasta(tempf,2,1);
       alnpro=translatealn(alndna);
      % viewseq(aln2)
       viewseqt(alndna.seq)
    end
end