function varargout = MBEGUI(varargin)
%MBEGUI - MBEToolbox GUI
%      MBEGUI, by itself, creates a new MBEGUI or raises the existing
%      singleton*.
%
%      H = MBEGUI returns the handle to a new MBEGUI or the handle to
%      the existing singleton*.
%
%      MBEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MBEGUI.M with the given input arguments.
%
%      MBEGUI('Property','Value',...) creates a new MBEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MBEGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MBEGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MBEGUI

% Last Modified by GUIDE v2.5 29-Dec-2012 20:05:57
global mbeversionstr
mbeversionstr='2.20';

% Molecular Biology and Evolution Toolbox (MBEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% Website: http://www.hku.hk/jamescai
% 
% $LastChangedDate: 2013-04-07 11:23:54 -0500 (Sun, 07 Apr 2013) $
% $LastChangedRevision: 524 $
% $LastChangedBy: jcai $



% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MBEGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MBEGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before MBEGUI is made visible.
function MBEGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MBEGUI (see VARARGIN)

% Choose default command line output for MBEGUI
global mbeguiolddir;
handles.output = hObject;

scrsz=get(0,'ScreenSize');
pos_act=get(gcf,'Position');
xr=scrsz(3)-pos_act(3);
xp=round(xr/2);
yr=scrsz(4)-pos_act(4);
yp=round(yr/2);
set(gcf,'position',[xp yp pos_act(3) pos_act(4)]);

if ~ispref('MBEToolbox','lastworkingdir')
   addpref('MBEToolbox','lastworkingdir',pwd)
else
   try
    cd(getpref('MBEToolbox','lastworkingdir',pwd));
   catch ME
    rmpref('MBEToolbox','lastworkingdir');
   end

end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MBEGUI wait for user response (see UIRESUME)
% uiwait(handles.MBEToolboxGUI);
% load 'MBEGUI.mat' MBEGUI_SaveDistance
%dirstr = fileparts(which(mfilename));
%cd(dirstr);

mbeguiolddir=pwd;
%cdmbe;

% deletetempfiles;
% clear global aln aln_ori
SetMenuStatus(handles);


% --- Outputs from this function are returned to the command line.
function varargout = MBEGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Sequences_Callback(hObject, eventdata, handles)
% hObject    handle to Sequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Distances_Callback(hObject, eventdata, handles)
% hObject    handle to Distances (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Phylogeny_Callback(hObject, eventdata, handles)
% hObject    handle to Phylogeny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Graph_Callback(hObject, eventdata, handles)
% hObject    handle to Graph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenPhylipFile_Callback(hObject, eventdata, handles)
% hObject    handle to OpenPhylipFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
	ButtonName=questdlg('What kind of Phylip format?', ...
			    'Select Phylip format', ...
			    'Interleaved','Sequential','Interleaved');
	switch ButtonName,
	    case 'Sequential',
           aln2=readphylip_s;
	    case 'Interleaved'
           aln2=readphylip_i;
	    otherwise
		return;
	end

if ~(isempty(aln2)),
    aln=aln2;
    SetMenuStatus(handles);
else
    return;
end

% --------------------------------------------------------------------
function OpenFASTAFile_Callback(hObject, eventdata, handles)
% hObject    handle to OpenFASTAFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
aln2=readfasta;
if ~(isempty(aln2)),
    aln=aln2;
    SetMenuStatus(handles);
else
    return;
end

% --------------------------------------------------------------------
function SaveFASTAFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFASTAFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
writefasta(aln);

% --------------------------------------------------------------------
function SavePhylipFile_Callback(hObject, eventdata, handles)
% hObject    handle to SavePhylipFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
writephylip_i(aln);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori mbeguiolddir mbeversionstr;
aln=[]; aln_ori=[];
cd(mbeguiolddir);
clear global aln aln_ori mbeversionstr
close;

% --------------------------------------------------------------------
function RemoveAln_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveAln (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
aln=[]; aln_ori=[];
clear global aln aln_ori
SetMenuStatus(handles);

% --------------------------------------------------------------------
function ViewSequences_Callback(hObject, eventdata, handles)
% hObject    handle to ViewSequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
	viewseq(aln);
catch ME
    errordlg(ME.message)
end



% --------------------------------------------------------------------
function ViewSequencesInfo_Callback(hObject, eventdata, handles)
% hObject    handle to ViewSequencesInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
aln
if (aln.seqtype==2),
    disp('*  Sequences in alignment (ALN) are codons (seqtype=2).')
    disp(sprintf('** Number of condons is %d.',size(aln.seq,2)./3));
end

% --------------------------------------------------------------------
function RemoveGapsInSequences_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveGapsInSequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
p=fileparts(which(mfilename));
cd(getpref('DesktopGenomics','workingdir',p));
cd('MBEToolbox');cd('addins');cd('trimAl'); cd('bin');
if (hasgap(aln))
	try
	    aln_ori=aln;
        writefasta(aln,'Input.fasta');
        trimcolumn;  % to trim the colunms with 'trimal' program %
        if isempty(aln.seq)
            disp('Empty Sequence.')
            aln
        else
            viewseq(aln);
        end

	    SetMenuStatus(handles);
	catch
	    errordlg(lasterr)
	end
else
    return;
end


% --------------------------------------------------------------------
function IncludeExcludeSequences_Callback(hObject, eventdata, handles)
% hObject    handle to IncludeExcludeSequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
[s,v] = choosebox('Name','Pick outgroup sequence(s)','PromptString',...
    'Sequences available:','SelectString','Selected sequences:',...
    'ListString',aln.seqnames);
if (v==1)
	try
	    aln_ori=aln;
	    aln.seq = aln.seq(s,:);
	    aln.seqnames=aln.seqnames(s);
	    SetMenuStatus(handles);
	catch
	    errordlg(lasterr)
	end
else
    return;
end;

% --------------------------------------------------------------------
function WorkOnAASequences_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnAASequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
	aln_ori=aln;
	[aln]=translatealn(aln);
	disp('Working on amino-acid sequences now!');
	set(handles.WorkOnAASequences,'enable','off');
	SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function WorkOnCodonPoistion1_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnCodonPoistion1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
	aln_ori=aln;
aln=extractpos(aln,1);
disp('Working on codon poistion 1 only now!');
SetMenuStatus(handles);

catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function WorkOnCodonPoistion2_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnCodonPoistion2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
aln_ori=aln;
aln=extractpos(aln,2);
disp('Working on codon poistion 2 only now!');
SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function WorkOnCodonPoistion3_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnCodonPoistion3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
aln_ori=aln;
aln=extractpos(aln,3);
disp('Working on codon poistion 3 only now!');
SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function WorkOnCodonPoistion12_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnCodonPoistion12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
aln_ori=aln;
aln=extractpos(aln,1,2);
disp('Working on codon poistion 1 and 2 only now!');
SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function WorkOnReverseStrand_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnReverseStrand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
aln_ori=aln;
aln=revseq(aln);
disp('Working on revserse of sequences now!');
SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function WorkOnRevcomStrand_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnRevcomStrand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
aln_ori=aln;
aln=revcomseq(aln);
disp('Working on reversed complement of sequences now!');
SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function WorkOnPloymorphicSites_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnPloymorphicSites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
aln_ori=aln;
aln=extractsegregatingsites(aln);
disp('Working on ploymorphic (segregating) sites only!');
SetMenuStatus(handles);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function WorkOnParsimonyInformativeSites_Callback(hObject, eventdata, handles)
% hObject    handle to WorkOnParsimonyInformativeSites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
aln_ori=aln;
[aln,informativesites]=extractinformativesites(aln);
disp('Working on Parsimony informative sites only!');
SetMenuStatus(handles);

catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function RestoreSequences_Callback(hObject, eventdata, handles)
% hObject    handle to RestoreSequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
try
    aln=aln_ori;
disp('Sequences restored!');
aln_ori=[];
SetMenuStatus(handles);

catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function ShowGCContent_Callback(hObject, eventdata, handles)
% hObject    handle to ShowGCContent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    [n,m]=size(aln.seq);
disp ('==========')
disp ('GC CONTENT')
disp ('==========')
for (k=1:n),
    [N,Nf,NN,NNf] = ntcomposition(aln.seq(k,:));
    GC=(N(2)+N(3))/(sum(N));
    fprintf('%s - %1.2f\n', aln.seqnames{k}, GC);
end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function ShowGC4Content_Callback(hObject, eventdata, handles)
% hObject    handle to ShowGC4Content (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    [V,P]=gc4(aln);
printmatrix(P,aln,2);
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function ShowNucleotideComposition_Callback(hObject, eventdata, handles)
% hObject    handle to ShowNucleotideComposition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    [n,m]=size(aln.seq);
disp ('======================')
disp ('Nucleotide Composition')
disp ('======================')
for (k=1:n),
    [N,Nf,NN,NNf] = ntcomposition(aln.seq(k,:));
    fprintf('%s\n', aln.seqnames{k});
    fprintf('A - %5d (%1.2f)\n', N(1),Nf(1));
    fprintf('C - %5d (%1.2f)\n', N(2),Nf(2));
    fprintf('G - %5d (%1.2f)\n', N(3),Nf(3));
    fprintf('T - %5d (%1.2f)\n\n', N(4),Nf(4));
end

catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function ShowCodonUsage_Callback(hObject, eventdata, handles)
% hObject    handle to ShowCodonUsage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    [V,P] = codonusage(aln,1);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function ShowCAI_Callback(hObject, eventdata, handles)
% hObject    handle to ShowCAI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    [CAI] = cai(aln,1);
printmatrix(CAI,aln);

catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function CalculateCodonVolatility_Callback(hObject, eventdata, handles)
% hObject    handle to calculateCodonVolatility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
disp ('================')
disp ('CODON VOLATILITY')
disp ('================')
[v,V] = codonvolatility(aln);
printmatrix(v,aln,5);
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function UPGMA_Callback(hObject, eventdata, handles)
% hObject    handle to UPGMA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
[D]=chooseDistance(aln);
if (isempty(D))
      return;
else
try
    figure;
    [topology] = plotupgma(D,char(aln.seqnames));
    %fprintf('remember: (although rooted by outgroup) this is an unrooted tree!');
catch ME
    errordlg(ME.message)
end

end


% --------------------------------------------------------------------
function NJ_Callback(hObject, eventdata, handles)
% hObject    handle to NJ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
[D]=chooseDistance(aln);
if (isempty(D))
      return;
else
try
    	figure;
	[anc,brnlen] = plotnjtree(D,1,char(aln.seqnames));
    fprintf('remember: this is an unrooted tree!');

catch ME
    errordlg(ME.message)
end

end

% --------------------------------------------------------------------
function DNAPars_Callback(hObject, eventdata, handles)
% hObject    handle to DNAPars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    runphylip(aln,'dnapars');
catch ME
    errordlg(ME.message)
end



% --------------------------------------------------------------------
function DNAML_Callback(hObject, eventdata, handles)
% hObject    handle to DNAML (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    runphylip(aln,'dnaml');
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function ProtPars_Callback(hObject, eventdata, handles)
% hObject    handle to ProtPars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    runphylip(aln,'protpars');
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function ProtML_Callback(hObject, eventdata, handles)
% hObject    handle to ProtML (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    runphylip(aln,'proml');
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function NTCompositionBarChart_Callback(hObject, eventdata, handles)
% hObject    handle to NTCompositionBarChart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    figure;
plotntcomposition(aln,2);
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function CorrespondenceAnalysisOfCodon_Callback(hObject, eventdata, handles)
% hObject    handle to CorrespondenceAnalysisOfCodon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    figure;
plotcorresp(aln,'rscu');
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function PlotDistanceMatrix_Callback(hObject, eventdata, handles)
% hObject    handle to PlotDistanceMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    D=dn_jc(aln);
	figure;
	%image(D.*100);
	pcolor(D)
	colormap(gray)
	%colormap(gray(2))
	axis ij
	axis square
	colorbar
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function PlotKaVsKs_Callback(hObject, eventdata, handles)
% hObject    handle to plotkavsks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    figure;
plotkavsks(aln);
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function PlotDistanceVsTransitionsTransversions_Callback(hObject, eventdata, handles)
% hObject    handle to PlotDistanceVsTransitionsTransversions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    figure;
plotdistvstrans(aln);

catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Plot3DZCurve_Callback(hObject, eventdata, handles)
% hObject    handle to Plot3DZCurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    figure;
plotzcurve(aln,'3d');

catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function PlotCumulativeGCProfile_Callback(hObject, eventdata, handles)
% hObject    handle to PlotCumulativeGCProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    figure;
plotzcurve(aln,'gc');

catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function PlotGCContent_Callback(hObject, eventdata, handles)
% hObject    handle to plotGCContent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    figure;
plotslidingwin(aln,'GC');
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function PlotAGContent_Callback(hObject, eventdata, handles)
% hObject    handle to plotAGContent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    figure;
plotslidingwin(aln,'AG');

catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function PlotGCDeviation_Callback(hObject, eventdata, handles)
% hObject    handle to plotGCDeviation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    figure;
plotslidingwin(aln,'GCDeviation');
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function PlotATDeviation_Callback(hObject, eventdata, handles)
% hObject    handle to plotATDeviation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
    figure;
plotslidingwin(aln,'ATDeviation');
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function HelpContents_Callback(hObject, eventdata, handles)
% hObject    handle to HelpContents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpwin mbetoolbox;

% --------------------------------------------------------------------
function NTDifference_Callback(hObject, eventdata, handles)
% hObject    handle to NTDifference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------xxx
function NTDifferenceD_Callback(hObject, eventdata, handles)
% hObject    handle to NTDifferenceD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[D,V,GAP]=dn_ntdiff(aln);
printmatrix(D,aln,0);
if i_WantToExportFile(handles),
    exportdismatrix(D,aln);
end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function NTDifferenceS_Callback(hObject, eventdata, handles)
% hObject    handle to NTDifferenceS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[P,Q]=countseqpq(aln);
printmatrix(P,aln,0);
if i_WantToExportFile(handles), exportdismatrix(P,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function NTDifferenceV_Callback(hObject, eventdata, handles)
% hObject    handle to NTDifferenceV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[P,Q]=countseqpq(aln);
printmatrix(Q,aln,0);
if i_WantToExportFile(handles), exportdismatrix(Q,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function NTDifferenceR_Callback(hObject, eventdata, handles)
% hObject    handle to NTDifferenceR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[P,Q]=countseqpq(aln);
D=P./Q;
printmatrix(D,aln,0);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function NTPDistance_Callback(hObject, eventdata, handles)
% hObject    handle to NTPDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function NTPDistanceD_Callback(hObject, eventdata, handles)
% hObject    handle to NTPDistanceD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln

try
[D,V]=dn_pdist(aln);
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function NTPDistanceS_Callback(hObject, eventdata, handles)
% hObject    handle to NTPDistanceS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln

try
[P,Q]=countseqpq(aln);[n,m]=size(aln.seq);
D=P./m;
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function NTPDistanceV_Callback(hObject, eventdata, handles)
% hObject    handle to NTPDistanceV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln

try
[P,Q]=countseqpq(aln);[n,m]=size(aln.seq);
D=Q./m;
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function NTPDistanceR_Callback(hObject, eventdata, handles)
% hObject    handle to NTPDistanceR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[P,Q]=countseqpq(aln);
D=P./Q;
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function AADifference_Callback(hObject, eventdata, handles)
% hObject    handle to AADifference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[D,V]=dp_aadiff(aln);
printmatrix(D,aln,0);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function AAPDistance_Callback(hObject, eventdata, handles)
% hObject    handle to AAPDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[D,V]=dp_pdist(aln);
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function AAPoissonDistance_Callback(hObject, eventdata, handles)
% hObject    handle to AAPoissonDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[D,V]=dp_poisson(aln);
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function AAGammaDistance_Callback(hObject, eventdata, handles)
% hObject    handle to AAGammaDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[D,V]=dp_gamma(aln);
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function JC69_Callback(hObject, eventdata, handles)
% hObject    handle to JC69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
D=dn_jc(aln);
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function K80_Callback(hObject, eventdata, handles)
% hObject    handle to K80 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
D=dn_k2p(aln);
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function TajimaNei84_Callback(hObject, eventdata, handles)
% hObject    handle to TajimaNei84 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
D=dn_tajima_nei84(aln);
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Felsenstein84_Callback(hObject, eventdata, handles)
% hObject    handle to f84 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
D=dn_f84(aln);
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Tamura92_Callback(hObject, eventdata, handles)
% hObject    handle to Tamura92 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
D=dn_tamura92(aln);
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function LogDetDistance_Callback(hObject, eventdata, handles)
% hObject    handle to LogDetDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
D=dn_logdet(aln);
printmatrix(D,aln,2);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function NtFreqEuclideanDistance_Callback(hObject, eventdata, handles)
% hObject    handle to NtFreqEuclideanDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[D]=dn_ntfreqdist(aln);
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function NG86_Callback(hObject, eventdata, handles)
% hObject    handle to NG86 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function NG86dS_Callback(hObject, eventdata, handles)
% hObject    handle to NG86dS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,VdS,VdN,St,Nt,nst,nat]=dc_ng86(aln);
printmatrix(dS,aln);
if i_WantToExportFile(handles), exportdismatrix(dS,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function NG86dN_Callback(hObject, eventdata, handles)
% hObject    handle to NG86dN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,VdS,VdN,St,Nt,nst,nat]=dc_ng86(aln);
printmatrix(dN,aln)
if i_WantToExportFile(handles), exportdismatrix(dN,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function NG86dNdS_Callback(hObject, eventdata, handles)
% hObject    handle to NG86dNdS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,VdS,VdN,St,Nt,nst,nat]=dc_ng86(aln);
printmatrix(dN_dS,aln)
if i_WantToExportFile(handles), exportdismatrix(dN_dS,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function NG86VdS_Callback(hObject, eventdata, handles)
% hObject    handle to NG86VdS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,VdS,VdN,St,Nt,nst,nat]=dc_ng86(aln);
printmatrix(VdS,aln)
if i_WantToExportFile(handles), exportdismatrix(VdS,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function NG86VdN_Callback(hObject, eventdata, handles)
% hObject    handle to NG86VdN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,VdS,VdN,St,Nt,nst,nat]=dc_ng86(aln);
printmatrix(VdN,aln)
if i_WantToExportFile(handles), exportdismatrix(VdN,aln); end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function NG86S_Callback(hObject, eventdata, handles)
% hObject    handle to NG86S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,VdS,VdN,St,Nt,nst,nat]=dc_ng86(aln);
printmatrix(St,aln)
if i_WantToExportFile(handles), exportdismatrix(St,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function NG86N_Callback(hObject, eventdata, handles)
% hObject    handle to NG86N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,VdS,VdN,St,Nt,nst,nat]=dc_ng86(aln);
printmatrix(Nt,aln)
if i_WantToExportFile(handles), exportdismatrix(Nt,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Li85_Callback(hObject, eventdata, handles)
% hObject    handle to Li85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Li85dS_Callback(hObject, eventdata, handles)
% hObject    handle to Li85dS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li85(aln);
printmatrix(dS,aln)
if i_WantToExportFile(handles), exportdismatrix(dS,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Li85dN_Callback(hObject, eventdata, handles)
% hObject    handle to Li85dN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li85(aln);
printmatrix(dN,aln)
if i_WantToExportFile(handles), exportdismatrix(dN,aln); end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function Li85dNdS_Callback(hObject, eventdata, handles)
% hObject    handle to Li85dNdS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li85(aln);
printmatrix(dN_dS,aln)
if i_WantToExportFile(handles), exportdismatrix(dN_dS,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Li85D4_Callback(hObject, eventdata, handles)
% hObject    handle to Li85D4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li85(aln);
printmatrix(D4,aln)
if i_WantToExportFile(handles), exportdismatrix(D4,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Li85D0_Callback(hObject, eventdata, handles)
% hObject    handle to Li85D0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li85(aln);
printmatrix(D0,aln)
if i_WantToExportFile(handles), exportdismatrix(D0,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Li85L4_Callback(hObject, eventdata, handles)
% hObject    handle to Li85L4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li85(aln);
printmatrix(L4,aln,1)
if i_WantToExportFile(handles), exportdismatrix(L4,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Li85L0_Callback(hObject, eventdata, handles)
% hObject    handle to Li85L0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li85(aln);
printmatrix(L0,aln,1)
if i_WantToExportFile(handles), exportdismatrix(L0,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Li93_Callback(hObject, eventdata, handles)
% hObject    handle to Li93 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Li93dS_Callback(hObject, eventdata, handles)
% hObject    handle to Li93dS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li93(aln);
printmatrix(dS,aln)
if i_WantToExportFile(handles), exportdismatrix(dS,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Li93dN_Callback(hObject, eventdata, handles)
% hObject    handle to Li93dN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li93(aln);
printmatrix(dN,aln)
if i_WantToExportFile(handles), exportdismatrix(dN,aln); end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function Li93dNdS_Callback(hObject, eventdata, handles)
% hObject    handle to Li93dNdS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li93(aln);
printmatrix(dN_dS,aln)
if i_WantToExportFile(handles), exportdismatrix(dN_dS,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Li93D4_Callback(hObject, eventdata, handles)
% hObject    handle to Li93D4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li93(aln);
printmatrix(D4,aln)
if i_WantToExportFile(handles), exportdismatrix(D4,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Li93D0_Callback(hObject, eventdata, handles)
% hObject    handle to Li93D0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li93(aln);
printmatrix(D0,aln)
if i_WantToExportFile(handles), exportdismatrix(D0,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Li93L4_Callback(hObject, eventdata, handles)
% hObject    handle to Li93L4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li93(aln);
printmatrix(L4,aln,1)
if i_WantToExportFile(handles), exportdismatrix(L4,aln); end
catch ME
    errordlg(ME.message)
end

% --------------------------------------------------------------------
function Li93L0_Callback(hObject, eventdata, handles)
% hObject    handle to Li93L0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS,L0,L2,L4,D0,D2,D4]=dc_li93(aln);
printmatrix(L0,aln,1);
if i_WantToExportFile(handles), exportdismatrix(L0,aln); end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function GY94dSdN_Callback(hObject, eventdata, handles)
% hObject    handle to MLdSdN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[dS,dN,dN_dS]=dc_codeml(aln);
disp('')
disp('==========================')
disp('dS')
printmatrix(dS,aln);
disp('')
disp('--------------------------')
disp('dN')
printmatrix(dN,aln);
disp('==========================')
% disp('dS/dN')
% printmatrix(dN_dS,aln);
if i_WantToExportFile(handles),
	exportdismatrix(dS,aln);
	exportdismatrix(dN,aln);
end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function KarlinGenomicSignature_Callback(hObject, eventdata, handles)
% hObject    handle to KarlinGenomicSignature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
[n,m]=size(aln.seq);
disp ('========================')
disp ('Karlin Genomic Signature')
disp ('========================')
for (k=1:n),
    [K] = karlinsig(aln.seq(k,:));
    fprintf('%s\n', aln.seqnames{k});
    fprintf('A %1.4f %1.4f %1.4f %1.4f\n', K(1,1),K(1,2),K(1,3),K(1,4));
    fprintf('C %1.4f %1.4f %1.4f %1.4f\n', K(2,1),K(2,2),K(2,3),K(2,4));
    fprintf('G %1.4f %1.4f %1.4f %1.4f\n', K(3,1),K(3,2),K(3,3),K(3,4));
    fprintf('T %1.4f %1.4f %1.4f %1.4f\n\n', K(4,1),K(4,2),K(4,3),K(4,4));
end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function PlotCummulativeKaKsCurves_Callback(hObject, eventdata, handles)
% hObject    handle to PlotCummulativeKaKsCurves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
figure;
[dS_cum,dN_cum,dS,dN] = plotslidingwinkaks(aln,20,2);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function PlotSlidingWindowKaKs_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSlidingWindowKaKs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
figure;
[dS_cum,dN_cum,dS,dN] = plotslidingwinkaks(aln,20,1);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function PlotDiPloMo_Callback(hObject, eventdata, handles)
% hObject    handle to plotdiplomo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
	figure;
	plotdiplomo(aln);
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function WantToExportDMToFile_Callback(hObject, eventdata, handles)
% hObject    handle to WantToExportDMToFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(hObject,'checked'),'on')
    set(hObject,'checked','off')
else
    set(hObject,'checked','on')
    set(handles.WantToExportDMToWS,'checked','off')
    set(handles.WantToExportDMToClipboard,'checked','off')    
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function WantToExportDMToClipboard_Callback(hObject, eventdata, handles)
% hObject    handle to WantToExportDMToClipboard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(hObject,'checked'),'on')
    set(hObject,'checked','off')
else
    set(hObject,'checked','on')
    set(handles.WantToExportDMToFile,'checked','off')
    set(handles.WantToExportDMToWS,'checked','off')    
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function ClustalW_Callback(hObject, eventdata, handles)
% hObject    handle to ClustalW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ClustalWalignDNAFile_Callback(hObject, eventdata, handles)
% hObject    handle to ClustalWalignDNAFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
AlnTemp = alignseqfile(1,0)
if (isempty(AlnTemp))
      return;
else
	aln=AlnTemp;
	viewseq(aln);
	SetMenuStatus(handles);
end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function ClustalWalignProteinFile_Callback(hObject, eventdata, handles)
% hObject    handle to ClustalWalignProteinFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
AlnTemp = alignseqfile(3,1)
if (isempty(AlnTemp))
      return;
else
	aln=AlnTemp;
	viewseq(aln);
	SetMenuStatus(handles);
end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function ClustalWalignCodingDNAFile_Callback(hObject, eventdata, handles)
% hObject    handle to ClustalWalignCodingDNAFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
try
AlnTemp = alignseqfile(2,1);
if (isempty(AlnTemp))
      return;
else
	aln=AlnTemp;
	viewseq(aln);
	SetMenuStatus(handles);
end
catch ME
    errordlg(ME.message)
end


% --------------------------------------------------------------------
function Demos_Callback(hObject, eventdata, handles)
% hObject    handle to Demo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
demo mbetoolbox


% --------------------------------------------------------------------
function AAJTTDistance_Callback(hObject, eventdata, handles)
% hObject    handle to AAJTTDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
[D]=dp_jtt(aln);
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end


% --------------------------------------------------------------------
function AADayhoffDistance_Callback(hObject, eventdata, handles)
% hObject    handle to AADayhoffDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
[D]=dp_dayhoff(aln);
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end

% --------------------------------------------------------------------
function AAWAGDistance_Callback(hObject, eventdata, handles)
% hObject    handle to AAWAGDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
[D]=dp_wag(aln);
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end


% --------------------------------------------------------------------
function SaveMatFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMatFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
    [filename, pathname, filterindex] = uiputfile( ...
       {'*.mat', 'MAT-file (*.mat)';
        '*.*',  'All Files (*.*)'}, ...
        'Save as');
	if ~(filename), return; end
	filename=[pathname,filename];

	if (filterindex==1)
	if (isempty(find(filename=='.')))
		filename=[filename,'.mat'];
	end
	end
save(filename,'aln');


% --------------------------------------------------------------------
function LoadMatFile_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMatFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
    [filename, pathname] = uigetfile( ...
       {'*.mat', 'MAT-file (*.mat)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a Mat file');
	if ~(filename), return; end
	filename=[pathname,filename];
alninmat = load(filename);

if ~(isfield(alninmat,'aln'))
    errordlg ('Invalid MAT file for alignment input')
    return;
else
    aln=alninmat.aln;
    SetMenuStatus(handles);
end


% --------------------------------------------------------------------
function HelpWebSite_Callback(hObject, eventdata, handles)
% hObject    handle to HelpWebSite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

web('http://bioinformatics.org/mbetoolbox/','-browser')

% --------------------------------------------------------------------
function CompEquTest_Callback(hObject, eventdata, handles)
% hObject    handle to CompEquTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
if (max(compequtest(aln))<=0.05),
        warndlg('Violation of compositional equilibrium detected!')
end



% --------------------------------------------------------------------
function SaveShadedHTMLFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveShadedHTMLFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
if (aln.seqtype~=3)
answ=questdlg('Save as Codon format?',...
'Output format','Yes','No','Cancel','No');
	switch (answ)
	    case ('Cancel')
		 return;
	    case ('Yes')
		writeshadedhtmlcodon(aln);
	    otherwise
		writeshadedhtml(aln);
	end
else
	writeshadedhtml(aln);
end


% --------------------------------------------------------------------
function SavePAMLFile_Callback(hObject, eventdata, handles)
% hObject    handle to SavePAMLFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
writepaml(aln);

% --------------------------------------------------------------------
function DNAVariability_Callback(hObject, eventdata, handles)
% hObject    handle to DNAVariability (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
figure; dnavariability(aln.seq,1);

% --------------------------------------------------------------------
function MatrixCircle_Callback(hObject, eventdata, handles)
% hObject    handle to MatrixCircle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MatrixCircleJTT_Callback(hObject, eventdata, handles)
% hObject    handle to MatrixCircleJTT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
m=modeljtt;
AA = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};
matrixcircle(m.R, AA, 'JTT')

% --------------------------------------------------------------------
function MatrixCircleGY94_Callback(hObject, eventdata, handles)
% hObject    handle to MatrixCircleGY94 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
icode=1;
[TABLE,C] = codontable;
	codon={'AAA' 'AAC' 'AAG' 'AAT' 'ACA' 'ACC' 'ACG' 'ACT' 'AGA' 'AGC' 'AGG' 'AGT' 'ATA' 'ATC' 'ATG' 'ATT'...
	      'CAA' 'CAC' 'CAG' 'CAT' 'CCA' 'CCC' 'CCG' 'CCT' 'CGA' 'CGC' 'CGG' 'CGT' 'CTA' 'CTC' 'CTG' 'CTT'...
	      'GAA' 'GAC' 'GAG' 'GAT' 'GCA' 'GCC' 'GCG' 'GCT' 'GGA' 'GGC' 'GGG' 'GGT' 'GTA' 'GTC' 'GTG' 'GTT'...
	      'TAA' 'TAC' 'TAG' 'TAT' 'TCA' 'TCC' 'TCG' 'TCT' 'TGA' 'TGC' 'TGG' 'TGT' 'TTA' 'TTC' 'TTG' 'TTT'};
stops=find(TABLE(icode,:)=='*');
codon(stops)=[];
m=modelgy94(0.5,3);
figure;
matrixcircle(m.R, codon, 'GY94 - omega=0.5, kappa=3')
xticklabel_rotate([1:61],90,codon,'interpreter','none')

% --------------------------------------------------------------------
function PlotEnhancedSlidingwindowKaKs_Callback(hObject, eventdata, handles)
% hObject    handle to PlotEnhancedSlidingwindowKaKs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
figure;
winsize=20;

id1=1; id2=2;
name1=aln.seqnames(id1);
name2=aln.seqnames(id2);
seq1=aln.seq(id1,:);
seq2=aln.seq(id2,:);

aln2=aln;
aln2.seqnames={char(name1), char(name2)};
aln2.seq=[seq1; seq2];

subplot(2,1,1); plotslidingwinkaks(aln2,winsize,1);
title([char(name1), ' -- ' ,char(name2), ' (window size = ', num2str(winsize), ')'])
subplot(2,1,2); plotslidingwinkaks(aln2,winsize,3);


% --------------------------------------------------------------------
function CompositionalComplexity_Callback(hObject, eventdata, handles)
% hObject    handle to CompositionalComplexity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
[n,m]=size(aln.seq);
disp ('========================')
disp ('COMPOSITIONAL COMPLEXITY')
disp ('========================')
for (k=1:n),
    K=compcomp(aln.seq);
    fprintf('%s - %1.4f\n', aln.seqnames{k}, K(k,1));
end


% --------------------------------------------------------------------
function CheckLatestVersion_Callback(hObject, eventdata, handles)
% hObject    handle to CheckLatestVersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mbeversionstr

try
	latestver = urlread('http://bioinformatics.org/mbetoolbox/LatestVersionTag');
catch
	warndlg('Can''t read URL.');
	return;
end

if (strcmp(latestver,mbeversionstr)==1),
    helpdlg('You are running the latest version.','No update available')
else
	ButtonName=questdlg(sprintf(['Download latest version %s?'], latestver), ...
			    'Check Latest Version', ...
			    'Go to website','Cancel','Go to website');
	switch ButtonName,
	    case 'Go to website',
            web('http://bioinformatics.org/mbetoolbox/','-browser')
        otherwise
		return;
	end
end
%web(['http://bioinformatics.org/mbetoolbox/mbetoolbox/version.asp?ver=',mbeversionstr],'-browser')


%% --------------------------------------------------------------------
%function CiteSoftware_Callback(hObject, eventdata, handles)
% hObject    handle to CiteSoftware (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%	ButtonName=questdlg('Looking at the reference?', ...
%			    'Cite MBEToolbox in your research', ...
%			    'Go to website','Cancel','Go to website');
%	switch ButtonName,
%	    case 'Go to website',
%           web('http://www.biomedcentral.com/1471-2105/6/64','-browser')
%	    otherwise
%		return;
%	end

% --------------------------------------------------------------------
function Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Stephens85_Callback(hObject, eventdata, handles)
% hObject    handle to Stephens85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln;
try
	stephens85(aln);
catch
	s=lasterror;
	warndlg(s.message);
	return;
end

% --------------------------------------------------------------------
function PLATO_Callback(hObject, eventdata, handles)
% hObject    handle to PLATO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln;
   [tree] = gettreedlg(aln);
   if(isempty(tree)),
       warndlg('No input tree.')
       return;
   end
   if(length(tree)<10),
       warndlg('Not a valid input tree.')
       return;
   end
      answer=questdlg('Do you want to continue?', ...
			    'Input tree loaded', ...
			    'Continue','Cancel','Continue');

      tree
      switch (lower(answer))
       case 'cancel'
           helpdlg('Action cancelled.')
           return;
       case 'continue'
           disp('Using JC69 model.')
           model=modeljc
	   figure;
	   plato(aln,tree,model);
       end


%warndlg('Function is under development. Sorry.')


% --------------------------------------------------------------------
function GenerateRandomSequences_Callback(hObject, eventdata, handles)
% hObject    handle to GenerateRandomSequences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ExportToWorkspace_Callback(hObject, eventdata, handles)
% hObject    handle to ExportToWorkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
%{
if ~isempty(aln),
assignin('base','aln',aln);
%disp('Sequences have been assigned as the variable, aln, in workspace.')
helpdlg('ALN assigned.','Assign varibles into workspace')
else
warndlg('Sequences is empty.','Assign varibles into workspace')
end
%}
try
    export2wsdlg({'Save Alignment (ALN) to Variable Named:'},...
                 {'aln'},{aln},...
                 'Export to Workspace');
catch ME
  errordlg(ME)
end



% --------------------------------------------------------------------
function SequenceInputdlg_Callback(hObject, eventdata, handles)
% hObject    handle to SequenceTextEditor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
   prompt={'FASTA format:'};
   def={''};
   dlgTitle='Sequences Input Dialog';
   lineNo=6;
   AddOpts.Resize='off';
   AddOpts.WindowStyle='modal';
   AddOpts.Interpreter='none';
   answer=inputdlg(prompt,dlgTitle,lineNo,def,AddOpts);

 if ~(isempty(answer)),
      seq=cellstr(answer{1});
      n=length(seq);
 if (n>=2&&~isempty(seq(2)))
       tempf = tempfilename;
       [fid,Msg] = fopen(tempf,'wt');
       if fid == -1, error(Msg); end
       for (k=1:n),
           fprintf(fid,'%s\n',char(seq(k)));
       end
       fclose(fid);

      answer=questdlg('Do you want to align input sequences?', ...
			    'Read Input Sequences', ...
			    'Yes','No','Yes');

      switch (lower(answer))
       case 'no'
	   aln=readfasta(tempf);
       case 'yes'
          [seqtype, geneticcode]=selectSeqTypeAndGeneticCode;
	   if (isempty(seqtype)||isempty(geneticcode)), aln=[]; return; end
	   aln = alignseqfile(seqtype,geneticcode,tempf);
       end

   SetMenuStatus(handles);
   aln
else
    helpdlg('Sorry, not valid sequences.')
end
end


% --------------------------------------------------------------------
function ListAA_Callback(hObject, eventdata, handles)
% hObject    handle to ListAA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listaa;


% --------------------------------------------------------------------
function RemoteDatabase_Callback(hObject, eventdata, handles)
% hObject    handle to RemoteDatabase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function SetMenuStatus(handles)
global aln aln_ori

if ~ispc
    set(handles.WantToExportDMToFile,'enable','off');
    set(handles.WantToExportDMToClipboard,'enable','off');
end

proteinonly = [handles.AAGammaDistance,...
    handles.AAPoissonDistance,...
    handles.AAPDistance,...
    handles.AADifference,...
    handles.AAJTTDistance,...
    handles.AADayhoffDistance,...
    handles.AAWAGDistance,...
    handles.AALGDistance,...
    handles.ProtPars,...
    handles.ProtML];

nucleotideonly=[handles.UPGMA,...
	handles.NJ,...
	handles.DNAPars,...
	handles.DNAML,...
	handles.NTDifference,...
	handles.NTPDistance,...
	handles.JC69,...
	handles.K80,...
	handles.Tamura92,...
	handles.TajimaNei84,...
	handles.Felsenstein84,...
	handles.LogDetDistance,...
	handles.NtFreqEuclideanDistance,...
	handles.NG86,...
	handles.Li85,...
	handles.Li93,...
	handles.GY94,...
    handles.GY94m,...
    handles.CompEquTest];

cannotemptyaln = [
        handles.RemoveAln,...
        handles.ExportToWorkspace,...
        handles.SaveFASTAFile,...
        handles.SavePhylipFile,...
        handles.SaveMEGAFile,...
        handles.SavePAMLFile,...
        handles.SaveMatFile,...
      handles.SaveShadedHTMLFile,...
      handles.ViewSequences,...
      handles.ViewSequencesInfo,...
      handles.RemoveGapsInSequences,...
      handles.WorkOnCodonPoistion1,...
      handles.WorkOnCodonPoistion2,...
      handles.WorkOnCodonPoistion3,...
      handles.WorkOnCodonPoistion12,...
      handles.WorkOnPloymorphicSites,...
      handles.WorkOnParsimonyInformativeSites,...
      handles.WorkOnRevcomStrand,...
      handles.WorkOnAASequences,...
      handles.RearrangeSequencesBySimilarity,...
      handles.IncludeExcludeSequences,...
      handles.IncludeExcludeSites,...
      handles.RestoreSequences,...
      handles.ShowGCContent,...
      handles.ShowGC4Content,...
      handles.ShowNucleotideComposition,...
      handles.ShowCodonUsage,...
      handles.ShowCAI,...
      handles.KarlinGenomicSignature,...
      handles.CompositionalComplexity,...
      handles.CalculateCodonVolatility,...
      handles.NTDifference,...
      handles.NTPDistance,...
      handles.AADifference,...
      handles.AAPDistance,...
      handles.AAPoissonDistance,...
      handles.AAGammaDistance,...
      handles.AADayhoffDistance,...
      handles.AAJTTDistance,...
      handles.AAWAGDistance,...
      handles.AALGDistance,...
      handles.JC69,...
      handles.K80,...
      handles.TajimaNei84,...
      handles.Tamura92,...
      handles.Felsenstein84,...
      handles.LogDetDistance,...
      handles.NtFreqEuclideanDistance,...
      handles.NG86,...
      handles.Li85,...
      handles.Li93,...
      handles.GY94,...
      handles.GY94m,...
      handles.CompEquTest,...
      handles.UPGMA,...
      handles.NJ,...
      handles.PhylipNeighbor,...
      handles.DNAPars,...
      handles.DNAML,...
      handles.ProtPars,...
      handles.ProtML,...
      handles.LikelihoodOfTree,...
      handles.Stephens85,...
      handles.Reticulate,...
      handles.PLATO,...
      handles.SiteSpecificRateML,...
      handles.SiteSpecificRateEB,...
      handles.NTCompositionBarChart,...
      handles.CorrespondenceAnalysisOfCodon,...
      handles.PlotDistanceMatrix,...
      handles.PlotDiPloMo,...
      handles.PlotDistanceVsTransitionsTransversions,...
      handles.PlotKaVsKs,...
      handles.PlotSlidingWindowKaKs,...
      handles.PlotEnhancedSlidingwindowKaKs,...
      handles.PlotCummulativeKaKsCurves,...
      handles.Plot3DZCurve,...
      handles.PlotCumulativeGCProfile,...
      handles.PlotGCContent,...
      handles.PlotAGContent,...
      handles.PlotGCDeviation,...
      handles.plotATDeviation,...
      handles.DNAVariability,...
      handles.DNASimilarity];

  if isempty(aln)
    set(cannotemptyaln,'enable','off');
  else
    set(cannotemptyaln,'enable','on');
    switch (aln.seqtype)
    case (1)
       set(proteinonly,'enable','off');
       set(nucleotideonly,'enable','on');
    case (2)
       set(proteinonly,'enable','off');
       set(handles.WorkOnAASequences,'enable','on');
       set(nucleotideonly,'enable','on');
    case (3)
       set(proteinonly,'enable','on');
       set(handles.WorkOnAASequences,'enable','off');
       set(nucleotideonly,'enable','off');
    end
    if ~(isempty(aln_ori)),
        set(handles.RestoreSequences,'enable','on');
    else
        set(handles.RestoreSequences,'enable','off');
    end
  end


% Update handles structure
% guidata(hObject, handles);

% --------------------------------------------------------------------
function Citation_Callback(hObject, eventdata, handles)
% hObject    handle to Citation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp(' ')
disp('Cai JJ, Smith DK, Xia X, Yuen KY.')
disp('MBEToolbox: a MATLAB toolbox for sequence data analysis in molecular biology and evolution.')
disp('BMC Bioinformatics. 2005 Mar 22;6(1):64.')
disp(' ')
disp('Cai JJ, Smith DK, Xia X, Yuen KY.')
disp('MBEToolbox 2.0: An enhanced version of a MATLAB toolbox for Molecular Biology and Evolution.')
disp('Evolutionary Bioinformatics Online 2006:2 189-192')
disp(' ')


% --------------------------------------------------------------------
function Reticulate_Callback(hObject, eventdata, handles)
% hObject    handle to Reticulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln;
try
	figure;
	reticulate(aln);
catch
	s=lasterror;
	warndlg(s.message);
	return;
end


% --------------------------------------------------------------------
function SiteSpecificRateML_Callback(hObject, eventdata, handles)
% hObject    handle to SiteSpecificRateML (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SiteSpecificERate(1)

% --------------------------------------------------------------------
function SiteSpecificRateEB_Callback(hObject, eventdata, handles)
% hObject    handle to SiteSpecificRateEB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SiteSpecificERate(2)

% --------------------------------------------------------------------
function SiteSpecificERate(method)
   global aln;
   [tree] = gettreedlg(aln);
   if(isempty(tree)),
       warndlg('No input tree.')
       return;
   end
   if(length(tree)<10),
       warndlg('Not a valid input tree.')
       return;
   end
	answer=questdlg('Do you want to continue?', ...
			    'Input tree loaded', ...
			    'Continue','Cancel','Continue');

      switch (lower(answer))
       case 'cancel'
           helpdlg('Action cancelled.')
           return;
       case 'continue'

             if (aln.seqtype~=3),
               model=modeljc;
               disp('Let''s use JC model')
             else
               model=modeljtt;
               disp('Let''s use JTT model')
             end

           [aln] = arrseq2tre(aln,tree);
           if (method==1),
                [rateraw,ratenorm]=rateofsites_ml(aln,tree,model,1);
                methodid='(ML method)';
           elseif (method==2),
                [rateraw,ratenorm]=rateofsites_eb(aln,tree,model,16,1);
                methodid='(EB method)';
           end
           figure;
           bar(ratenorm)
           title (['Site-specific Evolutionary Rate ',methodid]);
           ylabel('Normalised Relative Evolutionary Rate');
           xlabel('Site Position');
           axis([0 length(ratenorm)+1 min(ratenorm)*1.1 max(ratenorm)*1.1]);
  end


% --------------------------------------------------------------------
function GY94dSdNDirect_Callback(hObject, eventdata, handles)
% hObject    handle to MLdSdNDirect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln;
[s,v] = choosebox('Name','Select 2 sequences','PromptString',...
    'Sequences available:','SelectString','Selected 2 sequences:','ListString',aln.seqnames);
if (v==1)
    if (length(s)==2),
    currentaln=aln;
    currentaln.seq = aln.seq(s,:);
    currentaln.seqnames=aln.seqnames(s);
    [dS,dN,dN_dS]=dc_gy94(currentaln,1,2);
    %[dS,dN]=dc_ml(aln,s(1),s(2));
    disp(' ')
    disp('===================================')
    disp('        dS        dN     dN/dS')
    fprintf('%10.4f%10.4f%10.4f\n', dS, dN, dN_dS);
    disp('====================================')
    disp(' ')
    else
    warndlg('This functin is slow, please only select 2 sequences.')
    end
%    if i_WantToExportFile(handles),
%    	exportdismatrix(dS,currentaln);
%    	exportdismatrix(dN,currentaln);
%    end
end;


% --------------------------------------------------------------------
function About_Callback(hObject, eventdata, handles)
% hObject    handle to About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mbeversionstr
info{1}='';
info{2}='Molecular Biology & Evolution Toolbox';
info{3}='';
info{4}=['Version ',mbeversionstr];
info{5}='';
info{6}='Author: James Cai';
info{7}='';
x=date;
info{8}=sprintf('Copyright %s All Rights Reserved',x(end-3:end));
info{9}='';
info{10}='Please reference this software when using as part of';
info{11}='research.';
info{12} ='';
helpdlg(info,'About');


% --------------------------------------------------------------------
function PhylipNeighbor_Callback(hObject, eventdata, handles)
% hObject    handle to PhylipNeighbor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
%[outgrno]=selectOutgroup(aln);
%if ~(outgrno), return; end
outgrno=1;
[tree] = mbe_neighbor(aln,outgrno);
disp(tree)

if ~(isempty(tree)),
    answer=questdlg('Do you want to plot tree?', ...
			    'OUTTREE generated', ...
			    'Continue','Cancel','Continue');

      switch (lower(answer))
       case 'cancel'
           return;
       case 'continue'
   		cmd = 'mbe_neighbor.m';
		dirstr = fileparts(which(cmd));

        oldpath=pwd;
        cdmbe;
        cd 'addins';
        cd 'njplot';
        if (ispc),
                outtree = [dirstr,'\outtree'];
                system(['copy ', '"', outtree '"']);
                system(['njplot.exe outtree']);
        else
                outtree = [dirstr,'/outtree'];
                system(['./njplot ', '"', outtree, '"']);
        end
        cd(oldpath);


    end
end


% --------------------------------------------------------------------
function RandomNonCodingDNA_Callback(hObject, eventdata, handles)
% hObject    handle to RandomNonCodingDNA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
[n,m,names]=GetNumAndLen;
seqtype=1;
if (n>0&&m>0)
   aln.seqtype = seqtype;
   aln.geneticcode = 0;
   aln.seqnames = names;
   aln.seq=seqgen(n,m,seqtype);
   SetMenuStatus(handles);
   aln
else
    warndlg('Error')
end


% --------------------------------------------------------------------
function RandomCodingDNA_Callback(hObject, eventdata, handles)
% hObject    handle to RandomCodingDNA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
[n,m,names]=GetNumAndLen;
seqtype=2;
if (n>0&&m>0)
   aln.seqtype = seqtype;
   aln.geneticcode = 1;
   aln.seqnames = names;
   aln.seq=seqgen(n,m,seqtype);
   SetMenuStatus(handles);
   aln
else
    warndlg('Error')
end


% --------------------------------------------------------------------
function RandomProtein_Callback(hObject, eventdata, handles)
% hObject    handle to RandomProtein (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
[n,m,names]=GetNumAndLen;
seqtype=3;
if (n>0&&m>3)
   aln.seqtype = seqtype;
   aln.geneticcode = 1;
   aln.seqnames = names;
   aln.seq=seqgen(n,m,seqtype);
   SetMenuStatus(handles);
   aln
else
    warndlg('Error')
end


% --------------------------------------------------------------------
function [n,m,names]=GetNumAndLen
   prompt={'Number of sequences:','Length of sequences:'};
   def={'10','300'};
   dlgTitle='Generate random sequences';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);

 if ~(isempty(answer)),
   n=str2num(answer{1});
   m=str2num(answer{2});
   if(n>0),
       names={};
       for (k=1:n), names{k}=['seq',sprintf('%d',k)]; end
   end
   else
       n=[]; m=[]; names={};
 end
%names=reshape(num2cell(num2str(k(:), 'seq%1d'), 2), size(k));
% http://www.mathworks.com/matlabcentral/files/10006/mtx2charcell.m


% --------------------------------------------------------------------
function GY94mdSdNDirect_Callback(hObject, eventdata, handles)
% hObject    handle to GY94m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln;
[s,v] = choosebox('Name','Select 2 sequences','PromptString',...
    'Sequences available:','SelectString','Selected 2 sequences:','ListString',aln.seqnames);
if (v==1)
    if (length(s)==2),
    currentaln=aln;
    currentaln.seq = aln.seq(s,:);
    currentaln.seqnames=aln.seqnames(s);
    [dS,dN,dN_dS]=dc_gy94m(currentaln,1,2);
    %[dS,dN]=dc_ml(aln,s(1),s(2));
    disp(' ')
    disp('===================================')
    disp('        dS        dN     dN/dS')
    fprintf('%10.4f%10.4f%10.4f\n', dS, dN, dN_dS);
    disp('====================================')
    disp(' ')
    else
    warndlg('This functin is slow, please only select 2 sequences.')
    end
%    if i_WantToExportFile(handles),
%    	exportdismatrix(dS,currentaln);
%    	exportdismatrix(dN,currentaln);
%    end
end;


% --------------------------------------------------------------------
function MatrixCircleGY94m_Callback(hObject, eventdata, handles)
% hObject    handle to MatrixCircleGY94m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icode=1;
[TABLE,C] = codontable;
	codon={'AAA' 'AAC' 'AAG' 'AAT' 'ACA' 'ACC' 'ACG' 'ACT' 'AGA' 'AGC' 'AGG' 'AGT' 'ATA' 'ATC' 'ATG' 'ATT'...
	      'CAA' 'CAC' 'CAG' 'CAT' 'CCA' 'CCC' 'CCG' 'CCT' 'CGA' 'CGC' 'CGG' 'CGT' 'CTA' 'CTC' 'CTG' 'CTT'...
	      'GAA' 'GAC' 'GAG' 'GAT' 'GCA' 'GCC' 'GCG' 'GCT' 'GGA' 'GGC' 'GGG' 'GGT' 'GTA' 'GTC' 'GTG' 'GTT'...
	      'TAA' 'TAC' 'TAG' 'TAT' 'TCA' 'TCC' 'TCG' 'TCT' 'TGA' 'TGC' 'TGG' 'TGT' 'TTA' 'TTC' 'TTG' 'TTT'};
stops=find(TABLE(icode,:)=='*');
codon(stops)=[];
m=modelgy94m(0.5,3,2);
figure; matrixcircle(m.R, codon, 'GY94m - omega=0.5, kappa1=3, kapp2=2')
xticklabel_rotate([1:61],90,codon,'interpreter','none')


% --------------------------------------------------------------------
function DNASimilarity_Callback(hObject, eventdata, handles)
% hObject    handle to DNASimilarity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
figure; plotsimilarity(aln);


% --------------------------------------------------------------------
function IncludeExcludeSites_Callback(hObject, eventdata, handles)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori
[n,m]=size(aln.seq);
   prompt={'From site:','to:'};
   def={'1',num2str(m)};
   dlgTitle='Selectd Region';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);

 if ~(isempty(answer)),
     x=str2num(answer{1});
     y=str2num(answer{2});
    aln_ori=aln;
      aln.seq=aln.seq(:,[x:y]);
      SetMenuStatus(handles);
    aln
  end



% --------------------------------------------------------------------
function TranslationTable_Callback(hObject, eventdata, handles)
% hObject    handle to TranslationTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
showtranstable;


% --------------------------------------------------------------------
function OpenCLUSTALWFile_Callback(hObject, eventdata, handles)
% hObject    handle to OpenCLUSTALWFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
aln2=readclustal;
if ~(isempty(aln2)),
    aln=aln2;
    SetMenuStatus(handles);
else
    return;
end


% --------------------------------------------------------------------
function SaveCLUSTALWFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveCLUSTALWFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SaveMEGAFile_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMEGAFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
writemega(aln);




% --------------------------------------------------------------------
function CheckCodingDNAFile_Callback(hObject, eventdata, handles)
% hObject    handle to CheckCodingDNAFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
      [good,aln,filename]=checkcdsfile;
      if isempty(aln), return; end

      if (good),
          helpdlg('A validate CDS file.')
      else
      answer=questdlg('Do you want to view file?', ...
			    'Checking CDS File', ...
			    'Yes','No','Yes');
      switch (lower(answer))
       case 'no'
           helpdlg('Action cancelled.')
       case 'yes'
           dispfile(filename)
      end
   end


% --------------------------------------------------------------------
function SequenceTranslatedlg_Callback(hObject, eventdata, handles)
% hObject    handle to SequenceTranslatedlg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
translatedlg;


% --------------------------------------------------------------------
function RearrangeSequencesBySimilarity_Callback(hObject, eventdata, handles)
% hObject    handle to RearrangeSequencesByDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln aln_ori

[s,v] = choosebox('Name','Pick reference sequence','PromptString',...
    'Sequences available:','SelectString','Selected sequence:',...
    'ListString',aln.seqnames);
if (v==1)
    if length(s)>1
        warning('Only the first sequence is used as reference')
    end
    aln_ori=aln;
    aln=orderbysim(aln_ori,s(1));
    disp('Sequences in alignement have been rearraged (p-distance).')
    fprintf('Reference sequence: %s\n',aln_ori.seqnames{s(1)})
    fprintf('Before\t--->\tAfter\n')
    for k=1:length(aln.seqnames)
        fprintf('%4d %s\t%s\n',k,aln_ori.seqnames{k},aln.seqnames{k})
    end
    SetMenuStatus(handles);
else
    return;
end;



% --------------------------------------------------------------------
function LikelihoodOfTree_Callback(hObject, eventdata, handles)
% hObject    handle to LikelihoodOfTree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global aln;
   [tree] = gettreedlg(aln);
   if(isempty(tree)),
       warndlg('No input tree.')
       return;
   end
   if(length(tree)<10),
       warndlg('Not a valid input tree.')
       return;
   end
   tree
           disp('Using JC69 model.')
           model=modeljc
           lnL=treelike(aln,tree,model);
           assignin('base','aln',aln);
           assignin('base','tree',tree);
           assignin('base','model',model);
           assignin('base','lnL',lnL);
           fprintf('lnL = %f\n',lnL);

function GY94_Callback(hObject, eventdata, handles)
% hObject    handle to LikelihoodOfTree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function GY94m_Callback(hObject, eventdata, handles)
% hObject    handle to LikelihoodOfTree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function AALGDistance_Callback(hObject, eventdata, handles)
% hObject    handle to AALGDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global aln
%[D]=dp_lg(aln);
[D]=protdist(aln,'lg');
printmatrix(D,aln);
if i_WantToCopy(handles), num2clip(D); end
if i_WantToExportFile(handles), exportdismatrix(D,aln); end
if i_WantToExportWS(handles), exportdis2ws(D); end

% --------------------------------------------------------------------
function WantToExportDMToWS_Callback(hObject, eventdata, handles)
% hObject    handle to WantToExportDMToWS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(hObject,'checked'),'on')
    set(hObject,'checked','off')
else
    set(hObject,'checked','on')
    set(handles.WantToExportDMToFile,'checked','off')
    set(handles.WantToExportDMToClipboard,'checked','off')
end
guidata(hObject, handles);


% --------------------------------------------------------------------
% handles    structure with handles and user data (see GUIDATA)
function [yes]=i_WantToExportFile(handles)
if strcmp(get(handles.WantToExportDMToFile,'checked'),'on')
    yes=true;
else
    yes=false;
end

% --------------------------------------------------------------------
% handles    structure with handles and user data (see GUIDATA)
function [yes]=i_WantToCopy(handles)
if strcmp(get(handles.WantToExportDMToClipboard,'checked'),'on')
    yes=true;
else
    yes=false;
end

% --------------------------------------------------------------------
% handles    structure with handles and user data (see GUIDATA)
function [yes]=i_WantToExportWS(handles)
if strcmp(get(handles.WantToExportDMToWS,'checked'),'on')
    yes=true;
else
    yes=false;
end


function exportdis2ws(D)
try
    export2wsdlg({'Save Distance Matrix (D) to Variable Named:'},...
                 {'D'},{D},...
                 'Export to Workspace');
catch ME
  errordlg(ME)
end