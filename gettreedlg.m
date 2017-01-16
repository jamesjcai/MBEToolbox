function [tree] = gettreedlg(aln)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

	tree=[];

	if(aln.seqtype==3),
		seqtype='Protein';
	else
		seqtype='DNA';
	end

	ButtonName=questdlg('Please select a source of tree', ...
			    'We need an input tree', ...
			    'User tree','NJ tree','NJ tree');

	switch ButtonName,
	    case 'User tree', 
	        [tree] = i_usertree;
		return;
	     case 'NJ tree',	        
		[tree] = mbe_neighbor(aln);
	     otherwise
		return;
	end


function [answer] = i_usertree
   prompt={'Tree string:'};
   dlgTitle='Input Tree for Site-specific Rate Estimation';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo);

