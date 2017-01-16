function varargout = trimcolumn(varargin)
% TRIMCOLUMN MATLAB code for trimcolumn.fig
%      TRIMCOLUMN, by itself, creates a new TRIMCOLUMN or raises the existing
%      singleton*.
%
%      H = TRIMCOLUMN returns the handle to a new TRIMCOLUMN or the handle to
%      the existing singleton*.
%
%      TRIMCOLUMN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRIMCOLUMN.M with the given input arguments.
%
%      TRIMCOLUMN('Property','Value',...) creates a new TRIMCOLUMN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trimcolumn_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trimcolumn_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trimcolumn

% Last Modified by GUIDE v2.5 07-Feb-2013 21:38:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trimcolumn_OpeningFcn, ...
                   'gui_OutputFcn',  @trimcolumn_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before trimcolumn is made visible.
function trimcolumn_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trimcolumn (see VARARGIN)

% Choose default command line output for trimcolumn
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trimcolumn wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = trimcolumn_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in style1.
function style1_Callback(hObject, eventdata, handles)
% hObject    handle to style1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
system('trimal -in Input.fasta -out output1 -gt 1');
close;

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;


% --- Executes on button press in style2.
function style2_Callback(hObject, eventdata, handles)
% hObject    handle to style2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
system('trimal -in Input.fasta -out output2 -gt 0.8 -st 0.001');
close;

% --- Executes on button press in style3.
function style3_Callback(hObject, eventdata, handles)
% hObject    handle to style3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
system('trimal -in Input.fasta -out output3 -gt 0.8 -st 0.001 -cons 60');
close;
