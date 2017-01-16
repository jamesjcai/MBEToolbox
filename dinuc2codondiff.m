function [s,n,rad]=dinuc2codondiff(dinuc1,dinuc2)
% USAGE Example:
%
% DINU = {'AA','AC','AG','AT',...
%         'CA','CC','CG','CT',...
%         'GA','GC','GG','GT',...
%         'TA','TC','TG','TT'};
%
% for (i=1:length(DINU)),
%       for (j=i:length(DINU)),
% 	x=char(DINU{i});y=char(DINU{j});
% 	if (EditDist(x,y)==1)
% 	      fprintf(['%s <-> %s\n'],x,y);
% 	      dinuc2codondiff(x,y)
%             fprintf('\n');
% 	end
%       end
% end

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



c1=dinuc2codon(dinuc1);
c2=dinuc2codon(dinuc2);


[TABLE,CODON] = codontable;
icode=1;
T=TABLE(icode,:);

s=0; n=0;
rad=0;

printit=1;
switch (printit)
   case (1)
        fprintf(['\n%s <-> %s\n'],dinuc1,dinuc2);
	for (i=1:length(c1)),
	for (j=1:length(c2)),
		x=char(c1{i}); y=char(c2(j));
		a=findstr(x,dinuc1); b=findstr(y,dinuc2);
		if ((EditDist(x,y,2,2,1)==1)&&(sum(ismember(a,b))>0)),  % a single replacement allowed
			xx=codonise64(i_encode_n(x));
			yy=codonise64(i_encode_n(y));
			if ((T(xx)=='*')||(T(yy)=='*')),
				fprintf(['%s (%s) -- %s (%s), -\n'],c1{i},T(xx),c2{j},T(yy));
			else
				if (T(xx)==T(yy)),
					fprintf(['%s (%s) -- %s (%s), syn\n'],c1{i},T(xx),c2{j},T(yy));
					s=s+1;
				else
					cmark=' ';
					if (i_aaclass(T(xx))~=i_aaclass(T(yy))), cmark='^'; rad=rad+1; end
					fprintf(['%s (%s) -- %s (%s), nonsyn%s\n'],c1{i},T(xx),c2{j},T(yy),cmark);
					n=n+1;
				end
			end
		end
	end
	end
	if ((s+n)>0),
		fprintf(['Syn: %d (%2.2f); Nonsyn: %d (%2.2f)\n'],s,s/(s+n),n,n/(s+n));
		fprintf('Radical substitution: %d\n',rad);
	%fprintf(['%2.2f '],s/(s+n));
	end
     case (0)
	for (i=1:length(c1)),
	for (j=1:length(c2)),
		x=char(c1{i}); y=char(c2(j));
		a=findstr(x,dinuc1); b=findstr(y,dinuc2);
		if ((EditDist(x,y,2,2,1)==1)&&(sum(ismember(a,b))>0)),  % a single replacement allowed
			xx=codonise64(i_encode_n(x));
			yy=codonise64(i_encode_n(y));
			if ~((T(xx)=='*')||(T(yy)=='*')),
				if (T(xx)==T(yy)),
					s=s+1;
				else
					n=n+1;
				end
			end

		end
	end
	end

end



function [n] = i_aaclass(aa)
n=0;
classification1={{'R','H','K'};...    % positive
	         {'D','E'};...        % negative
      	         {'A','N','C','Q','G','I','L','M','F','P','S','T','W','Y','V'}};   % uncharged

classification2={{'D','E','H','K','R','N','Q','S','T'};...    % polar
	         {'L','I','V','M','F','Y','W','A','G','C','P'}};   % non-polar



classifior=classification1;
len=length(classifior);
for (k=1:len),
      if (ismember(aa,classifior{k}))
            n=k;
	    return
      end
end

% Classification by charge was made by dividing the
% amino acids into three categories: positive (R, H, K),
% negative (D, E), and uncharged (A, N, C, Q, G, I, L, M,
% F, P, S, T, W, Y, V).
%
% Classification by volume and polarity was made by
% dividing the amino acids into six categories: special (C),
% neutral and small (A, G, P, S, T), polar and relatively
% small (N, D, Q, E), polar and relatively large (R, H, K),
% nonpolar and relatively small (I, L, M, V), and nonpolar
% and relatively large (F, W, Y).