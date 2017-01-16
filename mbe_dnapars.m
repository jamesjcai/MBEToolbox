function varargout = mbe_dnapars(varargin)
%MBE_DNAPARS - M-file for mbe_dnapars.fig
%      MBE_DNAPARS, by itself, creates a new MBE_DNAPARS or raises the existing
%      singleton*.
%
%      H = MBE_DNAPARS returns the handle to a new MBE_DNAPARS or the handle to
%      the existing singleton*.
%
%      MBE_DNAPARS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MBE_DNAPARS.M with the given input arguments.
%
%      MBE_DNAPARS('Property','Value',...) creates a new MBE_DNAPARS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mbe_dnapars_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mbe_dnapars_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to RUN (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://bioinformatics.org/mbetoolbox/
% 
% $LastChangedDate: 2013-01-05 12:04:29 -0600 (Sat, 05 Jan 2013) $
% $LastChangedRevision: 327 $
% $LastChangedBy: jcai $



% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mbe_dnapars_OpeningFcn, ...
                   'gui_OutputFcn',  @mbe_dnapars_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if (nargin && ischar(varargin{1})),
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mbe_dnapars is made visible.
function mbe_dnapars_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mbe_dnapars (see VARARGIN)

% Choose default command line output for mbe_dnapars
handles.output = hObject;

scrsz=get(0,'ScreenSize');
pos_act=get(gcf,'Position');
xr=scrsz(3)-pos_act(3);
xp=round(xr/2);
yr=scrsz(4)-pos_act(4);
yp=round(yr/2);
set(gcf,'position',[xp yp pos_act(3) pos_act(4)]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mbe_dnapars wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% varargin{1};



% --- Outputs from this function are returned to the command line.
function varargout = mbe_dnapars_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Toption1.
function Toption1_Callback(hObject, eventdata, handles)
% hObject    handle to Toption1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Toption1
set(handles.Toption1, 'Value', 1);
set(handles.Toption2, 'Value', 0);
set(handles.Toption3, 'Value', 0);




% --- Executes on button press in Toption2.
function Toption2_Callback(hObject, eventdata, handles)
% hObject    handle to Toption2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Toption2
set(handles.Toption1, 'Value', 0);
set(handles.Toption2, 'Value', 1);
set(handles.Toption3, 'Value', 0);


% --- Executes on button press in Toption3.
function Toption3_Callback(hObject, eventdata, handles)
% hObject    handle to Toption3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Toption3
set(handles.Toption1, 'Value', 0);
set(handles.Toption2, 'Value', 0);
set(handles.Toption3, 'Value', 1);


% --- Executes on button press in RUN.
function RUN_Callback(hObject, eventdata, handles)
% hObject    handle to RUN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

optionstr = '';
i_Toption1 = get(handles.Toption1, 'Value');
i_Toption2 = get(handles.Toption2, 'Value');
i_Toption3 = get(handles.Toption3, 'Value');

if (i_Toption1==1)
    optionstr = [optionstr, ' -u1'];
elseif(i_Toption2==1)
    optionstr = [optionstr, ' -u2'];
elseif(i_Toption3==1)
    optionstr = [optionstr, ' -u3'];
end

s_maxtree = get(handles.maxtree, 'String');
    optionstr = [optionstr, ' -m', s_maxtree];

i_jumbleYES = get(handles.jumbleYES, 'Value');
i_jumbleNO = get(handles.jumbleNO, 'Value');
s_njumble = get(handles.njumble, 'String');

if (i_jumbleYES==1)
    optionstr = [optionstr, ' -j', s_njumble];
    optionstr = [optionstr, ' -s', num2str(i_oddseed)];
end

optionstr = [optionstr, ' -o', num2str(get(handles.outgrno,'Value'))];

i_threshYES = get(handles.threshYES, 'Value');
i_threshNO = get(handles.threshNO, 'Value');
s_threshold = get(handles.threshold, 'String');

if (i_threshYES==1)
    optionstr = [optionstr, ' -t', s_threshold];
end



i_transvpYES = get(handles.transvpYES, 'Value');
i_transvpNO = get(handles.transvpNO, 'Value');

if (i_transvpYES==1)
    optionstr = [optionstr, ' -v'];
end


if (exist('mbe_dnapars.fig','file')==2 && exist('mbe_dnapars.m','file')==2 )

dirstr=chdir2where('mbe_dnapars.m');
if (ispc)
		cmd = 'mbe_dnapars.exe';
		outtree = [dirstr,'\outtree'];
		outfile = [dirstr,'\outfile'];
else
		cmd = './mbe_dnapars';
		outtree = [dirstr,'/outtree'];
		outfile = [dirstr,'/outfile'];
end

% 	 if ((exist(outtree,'file')==2 && exist(outfile,'file')==2))
% 		option = questdlg('The file outfile/outtree already exist. Do you want to replace it?', ...
% 				   'File(s) Overwrite', ...
% 				  'Yes','No','No');
% 		switch option,
% 			case 'Yes',
%
% 			case 'No',
% 			return;
% 		end
% 	end



		cmd = [cmd, ' ', optionstr];
		disp(cmd)
		[s,w] = system(cmd);
if (s==0)
	disp(w)
	dispfile(outfile);
	answer=questdlg('Do you want to view output tree?','OUTTREE saved');

%dirstr = fileparts(which('MBEGUI'));
%if (ispc), sep='\'; else sep='/'; end
%njplotcmd = [dirstr,sep,'addins',sep,'njplot',sep,'njplot']
%njplotcmd=fullfile(dirstr,'addins','njplot','njplot');
olddir=pwd;
cdmbe;
cd('addins');
cd('njplot');
	switch (answer)
	    case ('Yes')
		     if (ispc)
			 system('njplot.exe outtree &');
		     else
			 system(',.njplot "', outtree, '" &');
		     end
	    otherwise

    end
cd(olddir);
else
    errordlg('Error occurred when running the third part program!');
    close;
end


	else
		errordlg('ERROR: File required not found.');
	end





% --- Executes on button press in CLOSE.
function CLOSE_Callback(hObject, eventdata, handles)
% hObject    handle to CLOSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;

% --- Executes during object creation, after setting all properties.
function maxtree_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxtree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function maxtree_Callback(hObject, eventdata, handles)
% hObject    handle to maxtree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxtree as text
%        str2double(get(hObject,'String')) returns contents of maxtree as a double


% --- Executes during object creation, after setting all properties.
function threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold as text
%        str2double(get(hObject,'String')) returns contents of threshold as a double


% --- Executes on button press in jumbleYES.
function jumbleYES_Callback(hObject, eventdata, handles)
% hObject    handle to jumbleYES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of jumbleYES
set(handles.jumbleYES, 'Value', 1);
set(handles.jumbleNO, 'Value', 0);
set(handles.njumble, 'Enable', 'On');


% --- Executes on button press in jumbleNO.
function jumbleNO_Callback(hObject, eventdata, handles)
% hObject    handle to jumbleNO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of jumbleNO
set(handles.jumbleYES, 'Value', 0);
set(handles.jumbleNO, 'Value', 1);
set(handles.njumble, 'Enable', 'Off');


% --- Executes during object creation, after setting all properties.
function outgrno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outgrno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in outgrno.
function outgrno_Callback(hObject, eventdata, handles)
% hObject    handle to outgrno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns outgrno contents as cell array
%        contents{get(hObject,'Value')} returns selected item from outgrno


% --- Executes on button press in transvpYES.
function transvpYES_Callback(hObject, eventdata, handles)
% hObject    handle to transvpYES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of transvpYES
set(handles.transvpYES, 'Value', 1);
set(handles.transvpNO, 'Value', 0);
set(handles.threshold, 'Enable', 'on');

% --- Executes on button press in transvpNO.
function transvpNO_Callback(hObject, eventdata, handles)
% hObject    handle to transvpNO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of transvpNO
set(handles.transvpYES, 'Value', 0);
set(handles.transvpNO, 'Value', 1);
set(handles.threshold, 'Enable', 'off');

% --- Executes on button press in threshYES.
function threshYES_Callback(hObject, eventdata, handles)
% hObject    handle to threshYES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of threshYES
set(handles.threshYES, 'Value', 1);
set(handles.threshNO, 'Value', 0);
set(handles.threshold, 'Enable','on');

% --- Executes on button press in threshNO.
function threshNO_Callback(hObject, eventdata, handles)
% hObject    handle to threshNO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of threshNO
set(handles.threshYES, 'Value', 0);
set(handles.threshNO, 'Value', 1);
set(handles.threshold, 'Enable','off');


% --- Executes during object creation, after setting all properties.
function njumble_CreateFcn(hObject, eventdata, handles)
% hObject    handle to njumble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function njumble_Callback(hObject, eventdata, handles)
% hObject    handle to njumble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of njumble as text
%        str2double(get(hObject,'String')) returns contents of njumble as a double


% --- Executes on button press in radiobutton14.
function radiobutton14_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton14


% --- Executes on button press in radiobutton15.
function radiobutton15_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton15



% --- Executes on button press in pushbuttonClose.
function pushbuttonClose_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;

