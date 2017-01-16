function varargout = mbe_proml(varargin)
%MBE_PROML - M-file for mbe_proml.fig
%      MBE_PROML, by itself, creates a new MBE_PROML or raises the existing
%      singleton*.
%
%      H = MBE_PROML returns the handle to a new MBE_PROML or the handle to
%      the existing singleton*.
%
%      MBE_PROML('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MBE_PROML.M with the given input arguments.
%
%      MBE_PROML('Property','Value',...) creates a new MBE_PROML or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mbe_proml_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mbe_proml_OpeningFcn via varargin.
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
                   'gui_OpeningFcn', @mbe_proml_OpeningFcn, ...
                   'gui_OutputFcn',  @mbe_proml_OutputFcn, ...
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


% --- Executes just before mbe_proml is made visible.
function mbe_proml_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mbe_proml (see VARARGIN)

% Choose default command line output for mbe_proml
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

% UIWAIT makes mbe_proml wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% varargin{1};



% --- Outputs from this function are returned to the command line.
function varargout = mbe_proml_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in RateOption1.
function RateOption1_Callback(hObject, eventdata, handles)
% hObject    handle to RateOption1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RateOption1
set(handles.RateOption1, 'Value', 1);
set(handles.RateOption2, 'Value', 0);
set(handles.RateOption3, 'Value', 0);
set(handles.alpha, 'Enable', 'Off');
set(handles.rcategs, 'Enable', 'Off');
set(handles.invarfrac, 'Enable', 'Off');

% --- Executes on button press in RateOption2.
function RateOption2_Callback(hObject, eventdata, handles)
% hObject    handle to RateOption2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RateOption2
set(handles.RateOption1, 'Value', 0);
set(handles.RateOption2, 'Value', 1);
set(handles.RateOption3, 'Value', 0);
set(handles.alpha, 'Enable', 'On');
set(handles.rcategs, 'Enable', 'On');
set(handles.invarfrac, 'Enable', 'off');

% --- Executes on button press in RateOption3.
function RateOption3_Callback(hObject, eventdata, handles)
% hObject    handle to RateOption3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RateOption3
set(handles.RateOption1, 'Value', 0);
set(handles.RateOption2, 'Value', 0);
set(handles.RateOption3, 'Value', 1);
set(handles.alpha, 'Enable', 'On');
set(handles.rcategs, 'Enable', 'On');
set(handles.invarfrac, 'Enable', 'On');


% --- Executes on button press in RUN.
function RUN_Callback(hObject, eventdata, handles)
% hObject    handle to RUN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

optionstr = '';


i_jttYES = get(handles.jttYES, 'Value');
i_jttNO = get(handles.jttNO, 'Value');
if (i_jttYES==1)
    optionstr = [optionstr, ' -m1'];
else
    optionstr = [optionstr, ' -m2'];
end


i_RateOption1 = get(handles.RateOption1, 'Value');
i_RateOption2 = get(handles.RateOption2, 'Value');
i_RateOption3 = get(handles.RateOption3, 'Value');
i_lambdaYES = get(handles.lambdaYES, 'Value');
i_lambdaNO = get(handles.lambdaNO, 'Value');
s_lambda = get(handles.lambda, 'String');

if (i_RateOption1==1)
    optionstr = [optionstr, ' -r1'];
elseif(i_RateOption2==1)
    optionstr = [optionstr, ' -r2'];
    optionstr = [optionstr, ' -p', get(handles.alpha, 'String')];
    optionstr = [optionstr, ' -c', num2str(get(handles.rcategs,'Value'))];
    if (i_lambdaYES==1)
        optionstr = [optionstr, ' -a', s_lambda];
    else
        optionstr = [optionstr, ' -a0'];
    end
elseif(i_RateOption3==1)
    optionstr = [optionstr, ' -r3'];
    optionstr = [optionstr, ' -p', get(handles.alpha, 'String')];
    optionstr = [optionstr, ' -c', num2str(get(handles.rcategs,'Value'))];
    optionstr = [optionstr, ' -i', get(handles.invarfrac, 'String')];
    if (i_lambdaYES==1)
        optionstr = [optionstr, ' -a', s_lambda];
    else
        optionstr = [optionstr, ' -a0'];
    end
end

i_jumbleYES = get(handles.jumbleYES, 'Value');
i_jumbleNO = get(handles.jumbleNO, 'Value');
s_njumble = get(handles.njumble, 'String');

if (i_jumbleYES==1)
    optionstr = [optionstr, ' -j', s_njumble];
    optionstr = [optionstr, ' -z', num2str(i_oddseed)];
else
    optionstr = [optionstr, ' -j0'];
end

optionstr = [optionstr, ' -o', num2str(get(handles.outgrno,'Value'))];

i_speedYES = get(handles.globalYES, 'Value');
i_speedNO = get(handles.globalNO, 'Value');
if (i_speedYES==1)
    optionstr = [optionstr, ' -s1'];
else
    optionstr = [optionstr, ' -s0'];
end

i_globalYES = get(handles.globalYES, 'Value');
i_globalNO = get(handles.globalNO, 'Value');
if (i_globalYES==1)
    optionstr = [optionstr, ' -g1'];
else
    optionstr = [optionstr, ' -g0'];
end




 if (exist('mbe_proml.fig','file')==2 && exist('mbe_proml.m','file')==2),

dirstr=chdir2where('mbe_proml.m');
if (ispc)
		cmd = 'mbe_proml.exe';
		outtree = [dirstr,'\outtree'];
		outfile = [dirstr,'\outfile'];
else
		cmd = '.\mbe_proml';
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
% 	   end
		cmd = [cmd, ' ', optionstr, ' '];
		disp(cmd)

[s,w] = system(cmd);
if (s==0)
	disp(w)
	dispfile(outfile);
	answer=questdlg('Do you want to view output tree?','OUTTREE saved');
	switch (answer)
	    case ('Yes')
		     if (ispc)
			 system(['njplot.exe ', '"', outtree, '"']);
		     else
			 system(['./njplot ', '"', outtree, '" &']);
		     end
	    otherwise

	end
else
    errordlg('Error occurred when running the third part program!');
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
function ttratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ttratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ttratio_Callback(hObject, eventdata, handles)
% hObject    handle to ttratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ttratio as text
%        str2double(get(hObject,'String')) returns contents of ttratio as a double


% --- Executes during object creation, after setting all properties.
function lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function lambda_Callback(hObject, eventdata, handles)
% hObject    handle to lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda as text
%        str2double(get(hObject,'String')) returns contents of lambda as a double


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


% --- Executes on button press in globalYES.
function globalYES_Callback(hObject, eventdata, handles)
% hObject    handle to globalYES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of globalYES
set(handles.globalYES, 'Value', 1);
set(handles.globalNO, 'Value', 0);
set(handles.lambda, 'Enable', 'on');

% --- Executes on button press in globalNO.
function globalNO_Callback(hObject, eventdata, handles)
% hObject    handle to globalNO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of globalNO
set(handles.globalYES, 'Value', 0);
set(handles.globalNO, 'Value', 1);
set(handles.lambda, 'Enable', 'off');

% --- Executes on button press in lambdaYES.
function lambdaYES_Callback(hObject, eventdata, handles)
% hObject    handle to lambdaYES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lambdaYES
set(handles.lambdaYES, 'Value', 1);
set(handles.lambdaNO, 'Value', 0);
set(handles.lambda, 'Enable','on');

% --- Executes on button press in lambdaNO.
function lambdaNO_Callback(hObject, eventdata, handles)
% hObject    handle to lambdaNO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lambdaNO
set(handles.lambdaYES, 'Value', 0);
set(handles.lambdaNO, 'Value', 1);
set(handles.lambda, 'Enable','off');


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


% --- Executes on button press in checkboxOUTFILE.
function checkboxOUTFILE_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxOUTFILE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxOUTFILE


% --- Executes on button press in checkboxTREEFILE.
function checkboxTREEFILE_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTREEFILE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxTREEFILE


% --- Executes on button press in pushbuttonView.
function pushbuttonView_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dirstr=chdir2where('mbe_proml.m');
if (ispc)
		outtree = [dirstr,'\outtree'];
		outfile = [dirstr,'\outfile'];
else
		outtree = [dirstr,'/outtree'];
		outfile = [dirstr,'/outfile'];
end

    i_checkboxOUTFILE = get(handles.checkboxOUTFILE, 'Value');
    i_checkboxTREEFILE = get(handles.checkboxTREEFILE, 'Value');

    if (i_checkboxOUTFILE==1), dispfile(outfile); end

     if (i_checkboxTREEFILE==1),
	     if (ispc)
		 system(['njplot.exe ', '"', outtree, '"']);
	     else
		 system(['./njplot ', '"', outtree, '" &']);
	     end
     end




% --- Executes on button press in pushbuttonClose.
function pushbuttonClose_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;


% --- Executes on button press in jttYES.
function jttYES_Callback(hObject, eventdata, handles)
% hObject    handle to jttYES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of jttYES
set(handles.jttYES, 'Value', 1);
set(handles.jttNO, 'Value', 0);


% --- Executes on button press in jttNO.
function jttNO_Callback(hObject, eventdata, handles)
% hObject    handle to jttNO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of jttNO
set(handles.jttYES, 'Value', 0);
set(handles.jttNO, 'Value', 1);



% --- Executes on button press in speedYES.
function speedYES_Callback(hObject, eventdata, handles)
% hObject    handle to speedYES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of speedYES
set(handles.speedYES, 'Value', 1);
set(handles.speedNO, 'Value', 0);


% --- Executes on button press in speedNO.
function speedNO_Callback(hObject, eventdata, handles)
% hObject    handle to speedNO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of speedNO
set(handles.speedYES, 'Value', 0);
set(handles.speedNO, 'Value', 1);



% --- Executes during object creation, after setting all properties.
function rcategs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rcategs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in rcategs.
function rcategs_Callback(hObject, eventdata, handles)
% hObject    handle to rcategs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns rcategs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from rcategs


% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function alpha_Callback(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha as text
%        str2double(get(hObject,'String')) returns contents of alpha as a double


% --- Executes during object creation, after setting all properties.
function invarfrac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to invarfrac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function invarfrac_Callback(hObject, eventdata, handles)
% hObject    handle to invarfrac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of invarfrac as text
%        str2double(get(hObject,'String')) returns contents of invarfrac as a double
