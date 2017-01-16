function [s_site,v_site,m_num,sn_num,sm_num]=countsegregatingsites2(aln)
%COUNTSEGREGATINGSITES2 - Counts segregating sites
%This function counts segregating sites in alignment(aln).
%
% Syntax: [s_site,v_site]=countsegregatingsites2(aln)
%
% Inputs:
%    aln   - Alignment structure
%
% Outputs:
%    s_site   - Number of segregating sites
%    v_site   - Number of valid sites = total length-gaps
%
% See also:

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $


if (isstruct(aln)), seq=aln.seq; else seq=aln; end
[n,m] = size(seq);
s_site = 0;
v_site = 0;
m_num = 0;    % number of mutation
sn_num = 0;   % number of singleton mutation
sm_num = 0;   % number of singleton mutation


for j=1:m
	thissite=seq(:,j);
	if (min(thissite)>0 & max(thissite) < 5),
		x=length(unique(thissite))-1;
		if (x>0),
			if (x==1 & i_issingleton(thissite)),
				sn_num=sn_num+1;
			end
			sm_num=sm_num+i_singleton_mut(thissite);
			m_num=m_num+x;
			s_site = s_site+1;
		end
	else
		v_site = v_site+1;	% take invalid char(>5 or <0) into account
	end
end
v_site = m-v_site;




function [y] = i_issingleton(site)
  y=0;
  [a,b,c]=unique(site);
  if (size(a,1)==2),
    if (sum(c==1)==1|sum(c==2)==1),
	y=1;
    end
  end


function [n] = i_singleton_mut(site)
  n=0;
  [a,b,c]=unique(site);
  p=size(a,1);
  for (k=1:p),
        n=n+(sum(c==k)==1);
  end