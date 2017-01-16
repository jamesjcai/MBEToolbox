function [tree] = gettreedlg2(aln)

	tree=[];

	if(aln.seqtype==3),
		seqtype='Protein';
	else
		seqtype='DNA';
	end

	ButtonName=questdlg('Do you want to input user tree?', ...
			    'Select source of input tree distance model', ...
			    'User tree',[seqtype,' Pars'],[seqtype, ' ML'],'User tree');

	switch ButtonName,
	    case 'User tree', 
	        [tree] = i_usertree;
		return;
	     case 'DNA Pars',	        
		cmd='mbetoolbox_dnapars.exe  -u1 -m100 -o1';
		%runaddin(aln,'dnapars',0)
	     case 'DNA ML',
		cmd='mbetoolbox_dnaml.exe  -t2.0 -f1 -r1 -j0 -o1 -s1 -g1';
		%runaddin(aln,'dnaml',0)
	     case 'Protein Pars',
		cmd='mbetoolbox_protpars.exe  -c1 -o1';
		%runaddin(aln,'protpars',0)
	     case 'Protein ML',
		cmd='mbetoolbox_proml.exe  -m1 -r1 -j0 -o1 -s1 -g1';
		%runaddin(aln,'proml',0)
	     otherwise
		return;
	end


	dirstr=chdir2where(['mbetoolbox_dnapars.exe']);
	writephylip_s(aln,[dirstr,'\infile']);

	[s,w] = system(cmd);
	if (s==0),
		disp(w)		
		tree=i_readouttree;
	end




function [answer] = i_usertree
   prompt={'Tree string:'};
   dlgTitle='Input Tree for Site-specific Rate Estimation';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo);


function [answer] = i_readouttree
	fid = fopen('outtree', 'r');
	answer = fscanf(fid, '%s');    % It has two rows now.
	fclose(fid);
	x=find(answer==';');
	if ~(isempty(x))
		answer=answer(1:x(1));
	end