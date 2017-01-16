function [res,tab,aln2] = paraduptab2(aln)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

a=rmcodongaps(aln)
c=a.seq;
s=codonise64(c);

[n,m]=size(s);

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

x=0;
res=zeros(1,m);
aasam=zeros(1,m);
for (k=1:m),
	st=s(:,k);
	st1=st(1); st2=st(2); st3=st(3); st4=st(4);
	if (st1==st2)			% same ref
	        x=x+1
		if (st3==st4)
			res(k)=1;	% same test
			aasam(k)=1;	% same test
		else
		      if (SYN(st3,st4))		% syn test
			    res(k)=2;
			    aasam(k)=1;
		      else			% nonsyn test
			    res(k)=3;
			    aasam(k)=2;
		      end
		end
	else
	      if (SYN(st1,st2))		% syn ref
			if (st3==st4)
				res(k)=4;	% same test
			        aasam(k)=1;
			else
			      if (SYN(st3,st4))		% syn test
				    res(k)=5;
         			    aasam(k)=1;
			      else			% nonsyn test
				    res(k)=6;
         			    aasam(k)=2;
			      end
			end
	      else			% nonsyn ref

    			if (st3==st4)
				res(k)=7;	% same test
        			aasam(k)=3;
			else
			      if (SYN(st3,st4))		% syn test
				    res(k)=8;
          			    aasam(k)=3;
			      else			% nonsyn test
				    res(k)=9;
        			aasam(k)=4;
			      end
			end

	      end
	end
%	res(k)
%	pause
	
end


tab=zeros(3);
tab(1,1)=sum(res==1);
tab(1,2)=sum(res==2);
tab(1,3)=sum(res==3);
tab(2,1)=sum(res==4);
tab(2,2)=sum(res==5);
tab(2,3)=sum(res==6);
tab(3,1)=sum(res==7);
tab(3,2)=sum(res==8);
tab(3,3)=sum(res==9);


disp('+=======================+')
disp('|      Same | Syn | NSyn|')
disp('+-----------------------+')
fprintf ('|Same |%5d|%5d|%5d|\n',tab(1,1),tab(1,2),tab(1,3));
disp('+-----------------------+')
fprintf ('|Syn  |%5d|%5d|%5d|\n',tab(2,1),tab(2,2),tab(2,3));
disp('+-----------------------+')
fprintf ('|NSyn |%5d|%5d|%5d|\n',tab(3,1),tab(3,2),tab(3,3));
disp('+=======================+')
disp('')
disp('')
disp('+===================+')
disp('|     | Same | Diff |')
disp('+-------------------+')
fprintf('|Same | %5d| %5d|\n',sum(aasam==1),sum(aasam==2));
disp('+-------------------+')
fprintf('|Diff | %5d| %5d|\n',sum(aasam==3),sum(aasam==4));
disp('+===================+')





cds=codonise64(c);

x1=cds(1,:)==cds(3,:);
x2=cds(2,:)~=cds(4,:);

x3=cds(1,:)~=cds(3,:);
x4=cds(2,:)==cds(4,:);

x12=x1.*x2;
x34=x3.*x4;
x1234=find((x12+x34)==1);

aln2=aln;
aln2.seq=decodonise64(cds(:,x1234));


