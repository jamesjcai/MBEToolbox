function [dS,dN,dN_dS]=dc_codeml(aln)
%DC_CODEML - Estimate dS and dN by using maximum likelihood algorithm
%This is a wrapper around the codeml program of PAML (Phylogenetic
%Analysis by Maximum Likelihood) package of Ziheng Yang.  See
%http://abacus.gene.ucl.ac.uk/software/paml.html for more information.
%
% Syntax: [dS,dN,dN_dS]=dc_codeml(aln)
%
% Inputs:
%    aln     - Alignment structure
%
% Outputs:
%    dS      - Synonymous substitution rate
%    dN      - Nonsynonymous substitution rate
%    dN_dS   - Ratio of dN and dS
%
% See also: DC_NEI_GOJOBORI86, DC_PAMILO_BIANCHI_LI93

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


warning off MATLAB:divideByZero

if (ispc)
	cmd = 'mbe_codeml.exe';
	sep='\';
else
	cmd = 'mbe_codeml';
	sep='/';
end
	dS=[];dN=[];dN_dS=[];


oldpath=pwd;
cdmbe;
cd 'addins';
cd 'paml';
if ~(exist(cmd,'file')==2) errordlg('ERROR: MBE_CODEML not found.'); return; end
%dirstr=chdir2where('dc_codeml.m');
%i_cleanfile;

writephylip_s(aln,'infile');

if ~(exist('mbe_codeml.ctl','file')==2)
	errordlg('ERROR: File MBE_CODEML.CTL required not found.');
	return;
else
%	h = waitbar(0.15,'Please wait...');
	if (ispc)
		cmd = [cmd, ' ', num2str(aln.geneticcode-1), ' '];
	else
		cmd = ['./',cmd, ' ', num2str(aln.geneticcode-1), ' '];
	end
		disp(['Running command: ',cmd])
		[s,w] = system(cmd);

	if s
	    error(w) % *** should add error processing
	   % close(h);
    else

%	        for i=15:95,
%    			waitbar(i/100,h)
%            end
%	        disp(w)
%    	    waitbar(1,h);
%            close(h);

	if (ispc),
		    dN2=dlmread('2ML.dN','');
	else
		    dN2=dlmread('./2ML.dN','');
	end

	    [m,n]=size(dN2);
	    if ~(m==n) dN2(:,end)=[]; end
	    dN=make_symmetric(dN2);

	if (ispc),
	    dS2=dlmread('2ML.dS','');
	else
	    dS2=dlmread('./2ML.dS','');
	end

	    [m,n]=size(dS2);	if ~(m==n) dS2(:,end)=[]; end
    	dS=make_symmetric(dS2);
	    dN_dS=diag2zeros(dN./dS);
	    dN_dS(isnan(dN_dS))=0;
	    dN_dS=diag2zeros(dN_dS);
	end
%	i_cleanfile;
end
cd(oldpath);




function [D] = make_symmetric(D2)
	[n,m]=size(D2);
	D=zeros(n,m);
	for (i=1:n),
	for (j=1:m),
		D(i,j)=D2(i,j);
		D(j,i)=D2(i,j);
	end
	end


function [D] = diag2zeros(D2)
	[n,m]=size(D2);
	D=D2;
	for (i=1:n),
	for (j=1:m),
		if (i==j)
			D(i,j)=0;
		end
	end
	end

function i_cleanfile()
	files={'2ML.dN','2ML.dS','2ML.t','2NG.dN','2NG.dS','2NG.t','infile','outfile','rst','rst1','rub'};
	for (k=1:length(files)),
	      file=char(files{k});
	      if (exist(file,'file')==2),
    		delete(file);
	      end
	end
