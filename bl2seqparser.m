function [p]=bl2seqparser(idx)

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $

    filename=sprintf('res/%d.res',idx);
    txt = textread(filename,'%s','delimiter','\n','whitespace','');
    % find first empty string in cell array, which occurs after the first
    % consensus line
    txt=strtrim(txt);
    mt = find(cellfun('isempty',txt));

% eliminate empty lines
txt(mt) = [];
qlen=txt{2};
qlen=str2num(qlen(2:strfind(qlen,' ')-1));
[id]=i_strfind(txt,'Query:');
if isempty(id)
    p=0;
    return
end

cline=zeros(1,qlen);
coord=[];
for k=1:length(id)
    %disp(txt{id(k)})
    a=txt{id(k)};
    [sid]=strfind(a,' ');
    s1=str2num(a(sid(1):sid(2)));
    s2=str2num(a(sid(end):end));
    coord=[coord;[s1 s2]];
    cline(s1:s2)=1;
end
%qlen
p=sum(cline)./qlen;






function [id]=i_strfind(txt,w)
id=[];
x=strfind(txt,w);
for k=1:length(x);
if ~isempty(x{k})
    id=[id,k];
end
end