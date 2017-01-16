function plottext(strText,xRatio,yRatio,color)
%Function plots with the Matlab function text(..) a string in the current
%figure. The text is placed at a position relative to the scale
%
%function plottext(strText,xRatio,yRatio,color)
%   strText:  Text, which should be plotted
%   xRatio:   x-position of the text relative to the x-scale [0 1]
%   yRatio;   y-position of the text relative to the y-scale [0 1]
%   color:    OPTIONAL, Color of the text, default is 'k' -> black
%
%EXAMPLE
%  The command
%  >> plottext('HELLO WORLD',0.5,0.5,'g')
%  plots the text 'HELLO WORLD' in the center of the current figure.
%  The color is green.

if nargin ==0
    help plottext
    return
elseif nargin<3
    help plottext
    return    
elseif nargin==3
    color = 'k';
end

x = get(gca,'xlim');
if strcmp(get(gca,'xscale'),'log')
    x = log10(x);
    dx = diff(x);
    x = 10^(x(1)+xRatio*dx);
else
    dx = diff(x);
    x = x(1)+xRatio*dx;    
end
y = get(gca,'ylim');
if strcmp(get(gca,'yscale'),'log')
    y = log10(y);
    dy = diff(y);
    y = 10^(y(1)+yRatio*dy);
else
    dy = diff(y);
    y = y(1)+yRatio*dy;    
end
h=text(x,y,strText,'fontsize',15);
set(h,'color',color)
