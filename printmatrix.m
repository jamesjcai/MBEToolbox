function printmatrix(X,aln,dec)
%PRINTMATRIX - Prints to the screen a printout of the matrix
%
% printmatrix(X,aln,dec)
% PRINTMATRIX prints to the screen a nice easy to read printout of the matrix X
% that can be copied and pasted into other applications (e.g. Excel).
%
% PRINTMATRIX(X); prints out the contents of X with 3 decimal places
%
% PRINTMATRIX(X,DEC); prints out the contents of X with DEC decimal places
%
%
% Written by Stephan W. Wegerich, SmartSignal Corp. August, 2001.
%
% if(nargin==1),dec=3;end

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


% [n,m]=size(aln.seqnames);
% for k=1:m
% fprintf(['%s\t'],char(aln.seqnames(n)));
% end
% fprintf('\n');

if(nargin<3),dec=3;end
if(nargin<2),aln.seqnames=strread(num2str(1:size(X,1)),'%s'); end

if(any(~isreal(X(:)))), error('Input Must be Real'); end

if(any(isnan(X(:))))
	disp('Raw output shown because NaN in result.');
	disp(' ');
	disp(X);
	return;
end

[N,M]=size(X);

ff = ceil(log10(max(max(abs(X)))))+dec+3;
fprintf('\n');

symmetry='nodiag';

	if ~(N==M)
		symmetry='yes';
	end
	
%symmetry='yes';
switch (symmetry)
    case ('yes')
	for i=1:N,
	    %name=char(aln.seqnames(i));
        name=aln.seqnames{i};
	    [x,len]=size(name);
		   if (len>10) 
			len=10;
			name = char(name(1:len));
		   elseif (len<10)
			name(len+1:10)=' ';	
		   end
	    fprintf(['%s '],name);


	    if (dec==0)
		fprintf(['%#',num2str(ff),'d '],X(i,:));
	    else
		fprintf(['%#',num2str(ff),'.',num2str(dec),'f '],X(i,:));
	    end	   
	    fprintf('\n');
	end
     case ('no')
	for i=1:N,
	    name=char(aln.seqnames(i));
	    [x,len]=size(name);
		   if (len>10)
			len=10;
			name = char(name(1:len));
		   elseif (len<10)
			name(len+1:10)=' ';
		   end
	    fprintf(['%s '],name);
	    if (i<=M)
		for j=1:i,
		    if (dec==0)
			fprintf(['%#',num2str(ff),'d '],X(i,j));
		    else
			fprintf(['%#',num2str(ff),'.',num2str(dec),'f '],X(i,j));
		    end
		end
	    end
        fprintf('\n');
	end

    case ('nodiag')
	for i=1:N,
	    name=char(aln.seqnames(i));
	    [x,len]=size(name);
		   if (len>10)
			len=10;
			name = char(name(1:len));
		   elseif (len<10)
			name(len+1:10)=' ';
		   end
	    fprintf(['%s '],name);

    if (i<=M)
	for j=1:i,
	if ~(i==j)
	    if (dec==0)
		fprintf(['%#',num2str(ff),'d '],X(i,j));
	    else
		fprintf(['%#',num2str(ff),'.',num2str(dec),'f '],X(i,j));
	    end
	end	
	end
      end
	    fprintf('\n');
    end
end