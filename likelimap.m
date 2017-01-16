function [prob,prob2]=likelimap(aln,model,flag)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

[n,m]=size(aln.seq);
if (n<4), error('LIKELIMAP requires at least 4 sequences.'); end
if (nargin<3), flag=1; end
if (nargin<2)
	if (aln.seqtype==3)
		model=modelwag;
	else
		model=modeljc;
	end
end
pick=nchoosek([1:n],4);
[pickn,dum]=size(pick);
prob=zeros(pickn,3);

fprintf('\n\nLIKELIHOOD MAPPING ANALYSIS\n\n');
fprintf('Number of quartets: %d  (all possible)\n', pickn);
fprintf('\nQuartet trees are based on approximate maximum likelihood values\n');
fprintf('using the selected model of substitution and rate heterogeneity.\n\n\n');

disp('Calculating likelihood...')
for (k=1:pickn),
    fprintf('Quartet %d of %d\n', k, pickn);
    s=aln.seq([pick(k,:)],:);
    prob(k,1)=treelike('((1,2),(3,4));',s,model);
    prob(k,2)=treelike('((1,3),(2,4));',s,model);
    prob(k,3)=treelike('((1,4),(2,3));',s,model);
end



	fprintf('\n\nLIKELIHOOD MAPPING STATISTICS\n\n');
	fprintf('Occupancies of the three areas 1, 2, 3:\n\n');
	fprintf('                        /\\\n');
	fprintf('                       /  \\\n');
	fprintf('                      /    \\\n');
	fprintf('                     /   1  \\\n');
	fprintf('                    / \\    / \\\n');
	fprintf('                   /   \\  /   \\\n');
	fprintf('                  /     \\/     \\\n');
	fprintf('                 /  3    :   2  \\\n');
	fprintf('                /        :       \\\n');
	fprintf('               /__________________\\\n');
	fprintf('\n');
	fprintf('Number of quartets in region 1: %d (= %.1f%%)\n', 5, 0.2);
	fprintf('Number of quartets in region 2: %d (= %.1f%%)\n', 5, 0.2);
	fprintf('Number of quartets in region 3: %d (= %.1f%%)\n\n', 5, 0.2);

prob2=i_qweight(prob);

if (flag),
   i_drawmap(prob2);
end




function i_drawmap(prob)
	ab=[bary2pt([1,0,0]);
	    bary2pt([0,1,0]);
	    bary2pt([0,0,1])];
	ab=ab';
	plot([ab(1,:),0], [ab(2,:),0])
	hold on

	[n,dum]=size(prob);
	for (k=1:n),
		pt=bary2pt(prob(k,:)./sum(prob(k,:)));
		plot(pt(1),pt(2),'r.')
	end



function [prob2] = i_qweight(prob)
	TP_MAX_EXP_DIFF=40;

	[n,m]=size(prob);
	prob2=zeros(n,m);

	for (k=1:n),
		p1=prob(k,1); p2=prob(k,2); p3=prob(k,3);
		%qweight = [p1,p2,p3];
		[qweight, qworder]=sort([-p1,-p2,-p3]);
		templog = qweight(2)-qweight(1);
			if(templog < -TP_MAX_EXP_DIFF),	% /* possible, since 1.0+exp(>36) == 1.0 */
			   qweight(2) = 0.0;
			else
			   qweight(2) = exp(templog);
			end

			templog = qweight(3)-qweight(2);
			if(templog < -TP_MAX_EXP_DIFF)	% /* possible, since 1.0+exp(>36) == 1.0 */
			   qweight(3)= 0.0;
			else
			   qweight(3) = exp(templog);
			end

			qweight(1) = 1.0;

			temp = qweight(1) + qweight(2) + qweight(3);
			prob2(k,1) = qweight(1)/temp;
			prob2(k,2) = qweight(2)/temp;
			prob2(k,3) = qweight(3)/temp;
	end
