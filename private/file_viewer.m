% This is a dialog box that pops up a little file viewer.
%-------------------------------------------------------------------------
% Bill Davidson, quellen@yahoo.com
% 16 Mar 2006
% It uses the GUI listbox control.
% It is stretchable.
% It ignores most control characters, so can view files containing binary.
% It is limited in file size to 500 lines (can be increased).
% It can work on '.mat' files, in which case it uses the diary, the
% command window, and a scratch file.
% It returns the figure handle.
%-------------------------------------------------------------------------
function varargout=file_viewer(varargin)

varargout{1}=NaN;
fpathfilename=varargin{1};
[pathstr,name,ext,versn]=fileparts(fpathfilename);
filename=[name ext];
usingscratch=0;
if strcmp(ext,'.mat')
    usingscratch=1;
    scratchfilename='foo.bar';
    s=load(fpathfilename);
    f=fieldnames(s);
    delete(scratchfilename);
    clc;
    format;
    diary(scratchfilename);
    for i=1:length(f)
        eval([char(f(i)) '=s.' char(f(i))])
    end
    diary off;
    clc;
    viewfilename=scratchfilename;
else
    viewfilename=fpathfilename;
end
if 2~=exist(viewfilename)
    return;
end

fig=openfig(mfilename,'reuse');
set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
handles=guihandles(fig);
guidata(fig,handles);

FigPos=get(0,'DefaultFigurePosition'); % where Matlab likes it to sit
OldUnits=get(fig, 'Units');
set(fig, 'Units', 'pixels');
OldPos=get(fig,'Position');
FigWidth=OldPos(3);
FigHeight=OldPos(4);
if isempty(gcbf)
    ScreenUnits=get(0,'Units');
    set(0,'Units','pixels');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);
    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits=get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos=get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2)=[(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
            (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
set(fig, 'Position', FigPos);
set(fig, 'Units', OldUnits);

maxlines=500; % max lines in file
text_lines=cell(1,maxlines);
fid=fopen(viewfilename,'r');
i=1;
filetoolarge=0;
while 1
    s=fgets(fid);
    if ~ischar(s), break, end
    s(s<32|s>127)=32; % ignore these characters
    text_lines(i)=cellstr(s);
    i=i+1;
    if i>maxlines
        filetoolarge=1;
        break;
    end
end
fclose(fid);
text_lines=text_lines(1:i-1);
handles.text_lines=text_lines;
guidata(handles.fooey,handles);
set(handles.file_viewer,'String',handles.text_lines,'Value',1);
if filetoolarge
    [s,err]=sprintf('File was truncated to %d lines',maxlines);warndlg(s,'file_viewer');
end
set(handles.fooey,'Name',filename);
varargout{1}=fig;
if usingscratch
    delete(scratchfilename);
    clc;
end
