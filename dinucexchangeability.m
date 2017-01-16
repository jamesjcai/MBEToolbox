function dinucexchangeability
%DINUCEXCHANGEABILITY - exchangeability of dinucleotides
%
%No reference, my idea.

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



DINU = {'AA','AC','AG','AT',...
        'CA','CC','CG','CT',...
        'GA','GC','GG','GT',...
        'TA','TC','TG','TT'};

data_syn=zeros(1,48);
data_nonsyn=zeros(1,48);
data_rad=zeros(1,48);
ct=0;

for (i=1:length(DINU)),
      for (j=i:length(DINU)),
	x=char(DINU{i});y=char(DINU{j});
	if (EditDist(x,y)==1)
	      ct=ct+1;
	      [d1,d2,d3]=dinuc2codondiff(x,y);
      	      data_syn(ct)=d1/(d1+d2);
	      data_nonsyn(ct)=d2/(d1+d2);
	      data_rad(ct)=d3;
	end
      end
end


%data_nonsynx=[1.00 0.57 1.00 1.00 1.00 1.00 1.00 0.50 1.00 1.00 1.00 1.00 0.71 1.00 1.00 1.00 1.00 1.00 0.50 0.50 0.50 1.00 1.00 0.50 0.50 1.00 1.00 0.50 1.00 1.00 1.00 0.75 0.71 0.57 0.71 1.00 0.75 0.50 1.00 0.75 1.00 1.00 0.50 0.50 0.50 0.71 0.50 0.71];
%data_synx=[0.00 0.43 0.00 0.00 0.00 0.00 0.00 0.50 0.00 0.00 0.00 0.00 0.29 0.00 0.00 0.00 0.00 0.00 0.50 0.50 0.50 0.00 0.00 0.50 0.50 0.00 0.00 0.50 0.00 0.00 0.00 0.25 0.29 0.43 0.29 0.00 0.25 0.50 0.00 0.25 0.00 0.00 0.50 0.50 0.50 0.29 0.50 0.29];

data1=data_syn;
data2=data_nonsyn;
data3=data_rad;

DINU = {'AA','AC','AG','AT',...
        'CA','CC','CG','CT',...
        'GA','GC','GG','GT',...
        'TA','TC','TG','TT'};
Map=[1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8; 1 9; 1 10; 1 11; 1 12; 1 13;...
     1 14; 1 15; 1 16; 2 3; 2 4; 2 5; 2 6; 2 7; 2 8; 2 9; 2 10; 2 11;...
     2 12; 2 13; 2 14; 2 15; 2 16; 3 4; 3 5; 3 6; 3 7; 3 8; 3 9; 3 10;...
     3 11; 3 12; 3 13; 3 14; 3 15; 3 16; 4 5; 4 6; 4 7; 4 8; 4 9; 4 10;...
     4 11; 4 12; 4 13; 4 14; 4 15; 4 16; 5 6; 5 7; 5 8; 5 9; 5 10; 5 11;...
     5 12; 5 13; 5 14; 5 15; 5 16; 6 7; 6 8; 6 9; 6 10; 6 11; 6 12; 6 13;...
     6 14; 6 15; 6 16; 7 8; 7 9; 7 10; 7 11; 7 12; 7 13; 7 14; 7 15; 7 16;...
     8 9; 8 10; 8 11; 8 12; 8 13; 8 14; 8 15; 8 16; 9 10; 9 11; 9 12;...
     9 13; 9 14; 9 15; 9 16; 10 11; 10 12; 10 13; 10 14; 10 15; 10 16;...
     11 12; 11 13; 11 14; 11 15; 11 16; 12 13; 12 14; 12 15; 12 16; 13 14;...
     13 15; 13 16; 14 15; 14 16; 15 16];

M1=zeros(16,16);
M2=zeros(16,16);
M3=ones(16,16)*-1;   % Radical label

dinuc1={'AA','AA','AA','AA','AA','AA','AC','AC','AC','AC','AC','AG',...
        'AG','AG','AG','AT','AT','AT','CA','CA','CA','CA','CA','CC',...
        'CC','CC','CC','CG','CG','CG','CT','CT','GA','GA','GA','GA',...
        'GC','GC','GC','GG','GG','GT','TA','TA','TA','TC','TC','TG'};
dinuc2={'AC','AG','AT','CA','GA','TA','AG','AT','CC','GC','TC','AT',...
        'CG','GG','TG','CT','GT','TT','CC','CG','CT','GA','TA','CG',...
        'CT','GC','TC','CT','GG','TG','GT','TT','GC','GG','GT','TA',...
        'GG','GT','TC','GT','TG','TT','TC','TG','TT','TG','TT','TT'};

for (k=1:48),
      x=dinuclise16(i_encode_n(char(dinuc1(k))));
      y=dinuclise16(i_encode_n(char(dinuc2(k))));
      M1(x,y)=data1(1,k);
      M2(x,y)=data2(1,k);
      M3(x,y)=data3(1,k);
end

tit='Dinucleotide Exchangeability';
i_matrixcircle_overlap(M1,M2,M3,DINU,tit)



function i_matrixcircle_overlap(M1,M2,M3,label,tit)

[n,m] = size(M1);
s1=sum(sum(triu(M1)));
s2=sum(sum(triu(M2)));
M1=(M1./(s1+s2)).*max(n,m)*2;
M2=(M2./(s1+s2)).*max(n,m)*2;

for (i=1:n),
for (j=i:n),
	if (i<j)
        nop=17;
        if (M1(i,j)==M2(i,j)), nop=10; end
              circle([i,j],M1(i,j),100,'r.');
	      circle([i,j],M2(i,j),nop,'.');  %nonsyn
	      if ~(ishold), hold on; end
	end
end
end

set(gca,'YTick',1:n)
set(gca,'YLim',[0 n+1])
set(gca,'YTickLabel',char(label));
%set(gca,'YDir','reverse')
set(gca,'XTick',1:n)
set(gca,'XLim',[0 n+1])
set(gca,'XTickLabel',char(label));
title(tit);
grid;
axis ij;             % OR set(gca,'YDir','reverse')
axis square;
xticklabel_rotate([1:16],90,label,'interpreter','none')
hold off;




figure;
for (i=1:n),
for (j=i:n),
	if (i<j)
		if (M3(i,j)>-1)
		      plottext(num2str(M3(i,j)),i,j,'b')
		      if ~(ishold), hold on; end
		end
	end
end
end

set(gca,'YTick',1:n)
set(gca,'YLim',[0 n+1])
set(gca,'YTickLabel',char(label));
%set(gca,'YDir','reverse')
set(gca,'XTick',1:n)
set(gca,'XLim',[0 n+1])
set(gca,'XTickLabel',char(label));
title('Radical substitution');
grid;
axis ij;             % OR set(gca,'YDir','reverse')
axis square;
xticklabel_rotate([1:16],90,label,'interpreter','none')
hold off;
