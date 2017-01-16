function [res,tab] = paraduptab(aln)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


a=rmcodongaps(aln)
s=a.seq;
c=codonise64(s);
[n,m]=size(c);

[TABLE,CODON] = codontable;
icode=1;
T=TABLE(icode,:);

SYN=zeros(64);
for (i=1:64),
for (j=1:64),
      if (T(i)==T(j))
	SYN(i,j)=1;
      end
end
end

res=zeros(1,m);
for (k=1:m),
	st=s(:,k);
	st1=st(1); st2=st(2); st3=st(3); st4=st(4);
	if (st1==st3)			% same ref
		if (st2==st4)
			res(k)=1;	% same test
		else
		      if (SYN(st2,st4))		% syn test
			    res(k)=2;
		      else			% nonsyn test
			    res(k)=3;
		      end
		end
	else
	      if (SYN(st1,st3))		% syn ref
			if (st2==st4)
				res(k)=4;	% same test
			else
			      if (SYN(st2,st4))		% syn test
				    res(k)=5;
			      else			% nonsyn test
				    res(k)=6;
			      end
			end
	      else			% nonsyn ref

    			if (st2==st4)
				res(k)=7;	% same test
			else
			      if (SYN(st2,st4))		% syn test
				    res(k)=8;
			      else			% nonsyn test
				    res(k)=9;
			      end
			end

	      end
	end
	
end


tab=zeros(3)
tab(1,1)=sum(res==1);
tab(1,2)=sum(res==2);
tab(1,3)=sum(res==3);
tab(2,1)=sum(res==4);
tab(2,2)=sum(res==5);
tab(2,3)=sum(res==6);
tab(3,1)=sum(res==7);
tab(3,2)=sum(res==8);
tab(3,3)=sum(res==9);
