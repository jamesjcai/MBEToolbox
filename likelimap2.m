function likelimap2(p1,p2,p3)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


ab=[bary2pt([1,0,0]);
    bary2pt([0,1,0]);
    bary2pt([0,0,1])];
ab=ab';
plot([ab(1,:),0], [ab(2,:),0])


a=[0 0.5 1]; b=[0 1 0];
P1=[a(1),b(1)];
P2=[a(2),b(2)];
P3=[a(3),b(3)];

%p1=1.3772; p2=0.4169; p3=5.2059;
pp=[p1 p2 p3]./sum([p1,p2,p3]);

PP=pp(1).*P1+pp(2).*P2+pp(3).*P3;



TP_MAX_EXP_DIFF=40;

%qweight = [p1,p2,p3];
[qweight, qworder]=sort([-p1,-p2,-p3])
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
	qweight(1) = qweight(1)/temp;
	qweight(2) = qweight(2)/temp;
	qweight(3) = qweight(3)/temp;

%	/* plot one point in likelihood mapping triangle */
	w1 = qweight(1);
	w2 = qweight(2);
	w3 = qweight(3);
	fprintf('%.10f tl %.10f tl dot\n', 0.5*w1 + w2, w1*0.8660254038);
	x1=0.5*w1 + w2;
	y1=w1*0.8660254038;


plot(a,b,PP(1),PP(2),'*',x1,y1,'rO');