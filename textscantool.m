function textscantool(varargin)
% Interactive tool to import and generate code to import data from delimited text
% files
%
% Example
%    textscantool

% $LastChangedDate: 2013-01-06 12:45:03 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 328 $
% $LastChangedBy: jcai $


% All possible format options for textscan plus string descriptions
formatOptions={...
    % 'formatchar', 'name', bytes'
    '%q',  'string', nan; % q for regular string
    '%f',  'double'  8;
    '%f32','single'  4;
    '%*q', 'ignore'  0;
    '%d8', 'int8'    1;
    '%d16','int16'   2;
    '%d32','int32'   4;
    '%d64','int64'   8;
    '%u8', 'uint8'   1;
    '%u16','uint16'  2;
    '%u32','uint32'  4;
    '%u64','uint64'  8;
    '%s',  'time'    8;... % Use s to designate time string
    };

% Default parameter values
numRowsTest=100;
numBlockRows=numRowsTest; % Number of rows to read in code gen when reading blocks
startingDataRow=1;
numHeaderLines=1; % Number of header lines
columnHeaderRowNumber=1; % Row for location of header columns
delimiter=','; % Delimiter

columnWidth=120; % Default width for columns in uitable
columnWidthLineStrings=1000;
numColumns=[]; % Preallocate number of columns (for use by other nested functions)
columnNames=cell(0); % Preallocate column names cell array of strings
columnNamesFormat=cell(0); % Preallocate strings
formatString=''; % Preallocate format sting to be used by textscan
bufferSize=4095; % Default buffer size to be used with textscan. Large for long rows
formatStringCell=cell(0); % Preallocate cell array of format strings
formatNameCell=cell(0); % Preallocate cell array of format string names
lineStrings=cell(0); % Preallocate cell array for lines of strings
pathName=pwd; % Path name
fid=[]; % File ID
dataCells2D=cell(0); % Preallocate 2D cells to load into uitable
codeLines=cell(0); % Prallocated cell array of generated lines of code
uitableJavaHandles=[]; % Handles to uitables
totalNumLinesText=[]; % Total number of lines in text file
dataStorageFormatOption=[]; % Data Storage option
importSize=[]; % Size in bytes
importTimePiece=[]; % Time to import
importTime=[]; % Time to import all estimate
importToBase=true; % Flag to tell import call back whether or not to inport to base
data=[]; % Store data between code gens to calculate size
GUIcolor=[0.9255    0.9137  0.8471];
selectedColumns=1;
numRowsCurrent=[]; % Number of data rows in file if less than 100 (small files)
generatedFunction='dataimport'; % Name of generated function

%% Initialze GUI layout
fig=openfig('textscantool.fig','new','invisible'); % Load GUI

properties=struct('menubar','none',...
    'NumberTitle','off',...
    'name','Textscan Tool: Select File, Specify Header and Delimiter Info',...
    'units','characters',...
    'closerequest',@mycloserequest);

set(fig,properties); % Set some properties

%% Previous, Next, Browse and Help buttons
openFileButtonHandle=findobj(fig,'Tag','openFileButton');
set(openFileButtonHandle,'callback',@fileopencallback);

helpPushButtonHandle=findobj(fig,'Tag','helpPushButton');
set(helpPushButtonHandle,'callback',@helpcallback);

pushButtonHandlesGUITags={'previousPushButton', 'nextPushButton'};
pushButtonHandlesCell=cellfun(@(x) findobj(fig,'Tag',x),pushButtonHandlesGUITags,'Uniform',false);
pushButtonHandles = cell2struct(pushButtonHandlesCell, pushButtonHandlesGUITags,2);

set(pushButtonHandles.previousPushButton,'callback',@previewpanel);
set(pushButtonHandles.previousPushButton,'Visible','on');
set(pushButtonHandles.nextPushButton,'callback',@formatpanel);
set(pushButtonHandles.nextPushButton,'Visible','on');

if nargin==1 % Remove
    fileName=varagin;
    set(pushButtonHandles.nextPushButton,'Enable','on');
else
    fileName='';
    set(pushButtonHandles.nextPushButton,'Enable','off');
end

%% Configure Preview panel uicontrols

% Find and name uicontrols
previewPanelGUITags={'fileInfoPanel','specifyPanel','previewPanel','fileNameText','fileSizeText',...
    'numHeaderLinesText','numHeaderLinesEditBox','columnHeadersCheckBox','columnHeadersRowText',...
    'columnHeadersRowEditBox','delimiterText','delimiterEditBox'};

previewPanelHandlesCell=cellfun(@(x) findobj(fig,'Tag',x),previewPanelGUITags,'Uniform',false);
previewPanelHandles = cell2struct(previewPanelHandlesCell, previewPanelGUITags,2);

% Configure uicontrols
set(previewPanelHandles.numHeaderLinesEditBox,'String',num2str(numHeaderLines));
set(previewPanelHandles.columnHeadersCheckBox,'Value',1);
set(previewPanelHandles.columnHeadersCheckBox,'callback',@checkboxcallback);
set(previewPanelHandles.columnHeadersRowEditBox,'String',num2str(columnHeaderRowNumber));
set(previewPanelHandles.delimiterEditBox,'String',delimiter);

% Disable initially (simplify with structfun as below?)
set(previewPanelHandles.numHeaderLinesEditBox,'Enable','off');
set(previewPanelHandles.columnHeadersCheckBox,'Enable','off');
set(previewPanelHandles.columnHeadersRowEditBox,'Enable','off');
set(previewPanelHandles.delimiterEditBox,'Enable','off');

% Make uitable
% *** Note: This old API for uitable is undocumented and unsupported as of MATLAB 7.5. The API may change ***
if verLessThan('MATLAB','7.6') % Before 8a
    uitableJavaHandles.lineStringsTable=uitable('Parent',fig,'Data',1, 'ColumnNames',{''},'gridcolor',GUIcolor);
else
    uitableJavaHandles.lineStringsTable=uitable('v0','Parent',fig,'Data',1, 'ColumnNames',{''},'gridcolor',GUIcolor);
end

uitableJavaHandles.lineStringsTable.setColumnWidth(columnWidthLineStrings);
set(uitableJavaHandles.lineStringsTable,'Units','characters')
set(uitableJavaHandles.lineStringsTable,'Editable',false);
pos=get(previewPanelHandles.previewPanel,'position'); % Position uitable relative to preview panel
distanceFromPreviewPanelY=-2.5;
distanceFromPreviewPanelX= 0.6;
uitableWidth=52;
uitableHeight=9;
set(uitableJavaHandles.lineStringsTable,'Position',[(pos(1)+distanceFromPreviewPanelX) (pos(2)+distanceFromPreviewPanelY) uitableWidth uitableHeight]);
set(uitableJavaHandles.lineStringsTable,'Visible',0); % Necessary?

% Set up call back to select row
rowHeaders=uitableJavaHandles.lineStringsTable.getTable.getParent.getParent.getRowHeader.getView;
set(rowHeaders, 'mouseclickedCallback', @rowselectcallback);
%rowHeaders.setRowHeader([]);

% Set invisible
structfun(@(x) set(x,'Visible','off'), previewPanelHandles);

%% Configure Format panel uicontrols

% File Name and size text
formatPanelGUITags={'dataTypePanel','dataViewPanel','timeFormatEdit'};
formatPanelHandlesCell=cellfun(@(x) findobj(fig,'Tag',x),formatPanelGUITags,'Uniform',false);
formatPanelHandles = cell2struct(formatPanelHandlesCell, formatPanelGUITags,2);

set(formatPanelHandles.dataTypePanel,'SelectionChangeFcn',@datatypepanelcallback)
set(formatPanelHandles.timeFormatEdit,'String','@(x) datenum(x, ''HH:MM:SS'')');  

% Make uitable
% *** Note: uitable is undocument and unsupported as of MATLAB 7.5. The API may change ***
if verLessThan('MATLAB','7.6') % Before 8a or later
    uitableJavaHandles.dataTable=uitable('Parent',fig,'Data',1, 'ColumnNames',{''});
else
    uitableJavaHandles.dataTable=uitable('v0','Parent',fig,'Data',1, 'ColumnNames',{''});
end
set(uitableJavaHandles.dataTable,'Units','characters')
pos=get(formatPanelHandles.dataViewPanel,'position'); % Position uitable relative to preview panel
distanceFromFormatPanelY=-2.5;
distanceFromFormatPanelX= 0.6;
set(uitableJavaHandles.dataTable,'Position',[(pos(1)+distanceFromFormatPanelX) (pos(2)+distanceFromFormatPanelY) (pos(3)-57) (pos(4)-9)]);
set(uitableJavaHandles.dataTable,'Editable',false);

% Set invisible
structfun(@(x) set(x,'Visible','off'), formatPanelHandles); % To do load invisible

%% Configure Import panel uicontrols

% Specify actions panel
importPanelGUITags={'importActionsPanel','howMuchImportPopup','dataStorageFormatPopup',...
    'generateCodeButton','importDataButton','dataStorageFormatText','howMuchText',...
    'numRowsText','numRowsEdit', 'countLinesButton', 'totalNumLinesText' 'importSizeText' 'inBlocksCheckBox',...
    'importTimeText','startingDataRowEdit'};

importPanelHandlesCell=cellfun(@(x) findobj(fig,'Tag',x),importPanelGUITags,'Uniform',false);
importPanelHandles = cell2struct(importPanelHandlesCell, importPanelGUITags,2);

% Configure
set(importPanelHandles.generateCodeButton,'callback',@generatefile);
set(importPanelHandles.importDataButton,'callback',@importDataCallBack);
set(importPanelHandles.dataStorageFormatPopup,'callback',@gencodelines);
set(importPanelHandles.inBlocksCheckBox,'callback',@gencodelines);
set(importPanelHandles.numRowsEdit,'String',num2str(numBlockRows));
set(importPanelHandles.numRowsEdit,'callback',@gencodelines);
set(importPanelHandles.startingDataRowEdit,'String',num2str(startingDataRow));
set(importPanelHandles.startingDataRowEdit,'callback',@gencodelines);
set(importPanelHandles.totalNumLinesText,'String','No. of lines in file: <unknown>');
set(importPanelHandles.countLinesButton,'callback',@countLinesCallBack);

% Set invisible
structfun(@(x) set(x,'Visible','off'), importPanelHandles);

%% Call preview panel
previewpanel()

%% End of Initialization

%% Preview panel
    function previewpanel(varargin) %#ok

        % Change title bar
        set(fig,'name','Textscan Tool: Select File');
        set(fig,'Visible','on');

        % On entering, hide/show appropriate panel uicontrols
        structfun(@(x) set(x,'Visible','on'), previewPanelHandles); % All not necessary?
        structfun(@(x) set(x,'Visible','off'), formatPanelHandles);

        if isfield(uitableJavaHandles,'dataTable') && ishandle(uitableJavaHandles.dataTable)
            datatablecallbackcontrol('off');
            set(uitableJavaHandles.dataTable,'Visible',0);
        end

        % Erase data
        dataCells2D=cell(0);

        % Configure Next button (is reset by format panel)
        set(pushButtonHandles.nextPushButton,'callback',@formatpanel);
        set(pushButtonHandles.nextPushButton,'String','Next >');

        % Make invisible previous push button
        set(pushButtonHandles.previousPushButton,'Visible','off'),

        if ~isempty(fileName) % if not defined
            populatepreviewpanel
        else
            set(previewPanelHandles.fileNameText,'String','');
            set(previewPanelHandles.fileSizeText,'String','Size: no file specified');
        end


    end

    function formatpanel(varargin) %#ok

        % Change title bar
        set(fig,'name','Textscan Tool: Specify Data Types of Columns or Choose Which to Ignore');

        % Make string panel uicontrols invisible
        set(uitableJavaHandles.lineStringsTable,'Visible',0); % Put with other handles below?

        % Hide/Show panel uicontrols
        structfun(@(x) set(x,'Visible','off'), previewPanelHandles);  % All not necessary?
        structfun(@(x) set(x,'Visible','on'), formatPanelHandles);
        structfun(@(x) set(x,'Visible','off'), importPanelHandles);

        % Configure Previous button
        set(pushButtonHandles.previousPushButton,'Visible','on'),
        set(pushButtonHandles.previousPushButton,'Enable','on');
        set(pushButtonHandles.previousPushButton,'callback',@previewpanel);
        set(pushButtonHandles.nextPushButton,'Visible','on');

        % Configure Next button
        set(pushButtonHandles.nextPushButton,'String','Next >');
        set(pushButtonHandles.nextPushButton,'callback',@importpanel);

        startingDataRow=1; % Reset inclase coming back from teh next panel
        
        if isempty(dataCells2D) % Text Files hasn't been read in formatted yet

            % Read values from uicontrols
            numHeaderLines=str2double(get(previewPanelHandles.numHeaderLinesEditBox,'String'));
            columnHeaderRowNumber=str2double(get(previewPanelHandles.columnHeadersRowEditBox,'String'));
            delimiter=get(previewPanelHandles.delimiterEditBox,'String');

            %% Read column headers if exists
            if get(previewPanelHandles.columnHeadersCheckBox,'Value'); % If has headers
                columnHeaderRow=lineStrings{columnHeaderRowNumber};
                columnNames=strread(columnHeaderRow,'%q','delimiter',delimiter)';
                numColumns=length(columnNames);

            else % If doesn't have headers
                columnHeaderRow=lineStrings{numHeaderLines+1}; % No headers
                columnNames=strread(columnHeaderRow,'%q','delimiter',delimiter)';
                numColumns=length(columnNames);
                extractSpaces=@(x) x(~(x==' ')); % Function to remove spaces from a string
                columnNameNumberStrings=cellstr(num2str((1:numColumns)'))'; % String numbers (may have spaces at start)
                columnNameNumberStrings=cellfun(extractSpaces,columnNameNumberStrings,'UniformOutput',false); % Remore spaces
                columnNames=cellfun(@(x) ['column' x],columnNameNumberStrings,'UniformOutput',false); % Make names that are just numbers
            end

            %% Format string Cell
            formatStringCell=repmat({'%q'},1,numColumns);
            formatNameCell=cellfun(@(x) formatChar2Name(x),formatStringCell,'UniformOutput',false);

            %% Data uitable

            % Read
            fseek(fid,0,'bof'); % rewind
            formatString=cell2mat(formatStringCell);
            dataCells=textscan(fid,formatString,numRowsTest,'delimiter',delimiter,'headerlines',numHeaderLines+(startingDataRow-1)); % Read strings delimited by a carriage return

            % Convert to cell of strings
            dataCells2D=textscandatato2Dcell(dataCells);

            % Make column names
            columnNamesFormat=cellfun(@(x,y) sprintf([x ' (' y ')']),...
                columnNames,formatNameCell,'UniformOutput',false);

            % Update uitable
            uitableJavaHandles.dataTable.setData(dataCells2D);
            uitableJavaHandles.dataTable.setColumnNames(columnNamesFormat);
            uitableJavaHandles.dataTable.setColumnWidth(columnWidth);
            set(uitableJavaHandles.dataTable,'Visible',1);

        else % Text file has already been read in formatted
            set(uitableJavaHandles.dataTable,'Visible',1);
        end

        datatablecallbackcontrol('on');
        dataTableCallback();

    end

    function importpanel(varargin) %#ok

        % Change title bar
        set(fig,'name','Textscan Tool: Specify How Code is to be Generated or Data Imported');

        % Is small file
        maxDataLines=totalNumLinesText-numHeaderLines;
        if (maxDataLines<100)
            numBlockRows=maxDataLines;
            set(importPanelHandles.numRowsEdit,'String',num2str(numBlockRows));
            set(importPanelHandles.totalNumLinesText,'String',['No. of lines in file:' int2str(totalNumLinesText)]);
        end

        % Make format panel uitable invisible
        datatablecallbackcontrol('off');
        set(uitableJavaHandles.dataTable,'Visible',0);

        % Hide/Show panel uicontrols
        structfun(@(x) set(x,'Visible','off'), formatPanelHandles);
        structfun(@(x) set(x,'Visible','on'), importPanelHandles);

        % Configure Previous  button
        set(pushButtonHandles.previousPushButton,'Enable','on');
        set(pushButtonHandles.previousPushButton,'callback',@formatpanel);

        % Configure Next button
        set(pushButtonHandles.nextPushButton,'Visible','off');

        % Gen code lines
        gencodelines
    end


    function checkboxcallback(varargin) %#ok

        if get(previewPanelHandles.columnHeadersCheckBox,'Value'); % If Column headers
            set(previewPanelHandles.columnHeadersRowEditBox,'Enable','on'); % Enabel edit box
        else
            set(previewPanelHandles.columnHeadersRowEditBox,'Enable','off'); % Otherwise disable box
        end
    end

    function datatypepanelcallback(buttonGroupHandle,changeInfo)

        % Column selection
        selectedColumns=1+uitableJavaHandles.dataTable.Table.getSelectedColumns;

        % Get new format char
        selectedHandle=changeInfo.NewValue;
        previousHandle=changeInfo.OldValue;
        fmt=get(selectedHandle,'Tag');
        formatName=fmt(1:end-5); % Remove radio string
        formatChar=formatName2Char(formatName); % Get \ format char

        % Make new format string cell array
        oldFormatString=formatString; % Take not of previous
        newFormatStringCell=formatStringCell;
        newFormatStringCell(selectedColumns')={formatChar};
        newformatString=cell2mat(newFormatStringCell);

        fseek(fid,0,'bof'); % Rewind

        % Read data
        try
            % Read data from file with requested format
            dataCell=textscan(fid,newformatString,numRowsTest,'delimiter',delimiter,'headerlines',numHeaderLines + (startingDataRow-1)); % Read strings delimited by a carriage return

            % Has succesfuly read in without errors at this point

            if isempty(dataCell) % If all ignores
                numRowsCurrent=numRowsTest; % Use default number of rows (to construct blanks)
            else
                numRowsCurrent=length(dataCell{1}); % Use current length
            end

            % If any time/date strings, convert to datenums
            found=strcmp(newFormatStringCell,'%s'); % Find time strings
            if any(found)
                stringFnc=get(formatPanelHandles.timeFormatEdit,'String');
                try
                    dateNumFncHandle=eval(stringFnc); % Evaluate function Handle
                    dataCell(found)=cellfun(dateNumFncHandle,dataCell(found),'UniformOutput',false);
                catch
                    errordlg('Invalid anonymous function for date time conversion')
                end
            end

            % If any 'ignores' in format string, replace with blanks
            found=strcmp(newFormatStringCell,'%*q'); % Find ignored colmuns
            if any(found)
                blankColumn={repmat({' '},numRowsCurrent,1)};
                dataCellNew(~found)=dataCell; % Copy in not ignored data
                dataCellNew(found)=blankColumn; % Pad in ignored columns with blanks
                dataCell=dataCellNew; % Replace data structure
            end

            % Convert 1D cell array to 2D cell array 
            dataCell2D=textscandatato2Dcell(dataCell);

            % Put data back into table
            uitableJavaHandles.dataTable.setData(dataCell2D);
            formatStringCell=newFormatStringCell;
            formatNameCell=cellfun(@(x) formatChar2Name(x),formatStringCell,'UniformOutput',false);

            % Update Column headers
            columnNamesFormat=cellfun(@(x,y) sprintf([x ' (' y ')']),columnNames,formatNameCell,'UniformOutput',false);

            datatablecallbackcontrol('off');
            set(uitableJavaHandles.dataTable,'ColumnNames',columnNamesFormat); %Put in
            datatablecallbackcontrol('on');

            % Reset column widths
            uitableJavaHandles.dataTable.setColumnWidth(columnWidth);

            % Set new format string
            formatString=newformatString;

        catch % If errors out
            if length(selectedColumns)==1
                warndlg(['The selected column ('  int2str(selectedColumns) ...
                    ') cannot be read as a ' formatName '. Try a different format.' ],...
                    'Invalid format');
            else
                warndlg(['Some of the selected columns (' int2str(selectedColumns') ...
                    ') cannot be read as ' formatName 's. Try a different format.' ],...
                    'Invalid format');
            end

            % Put button back
            set(buttonGroupHandle,'SelectedObject',previousHandle);
            formatString=oldFormatString;
        end

    end

    function charFormat=formatName2Char(nameFormat)
        % Converts 'single' to 'f32'
        found=strcmp(formatOptions(:,2), nameFormat);
        charFormat=formatOptions{found,1};
    end

    function nameFormat=formatChar2Name(charFormat)
        % Converts  'f32' to 'single'
        found=strcmp(formatOptions(:,1), charFormat);
        nameFormat=formatOptions{found,2};
    end

    function bytes=formatName2Bytes(formatName)
        % Converts  'f32' to 'single'
        found=strcmp(formatOptions(:,2), formatName);
        bytes=formatOptions{found,3};
    end

    function dataTableCallback(varargin) %#ok
        % When data table in format panel is selected, this function selects
        % the correct radio button if all columne are the same data type and
        % leaves it unselected if they are not.

        selectedColumns=1+uitableJavaHandles.dataTable.Table.getSelectedColumns; % Column selection
        currentFormatNames=formatNameCell(selectedColumns); % Get data type format names

        sameType=unique(currentFormatNames); % Find uniques data types
        if length(sameType)==1 % If all same
            rightButton=findobj(fig,'Tag',[sameType{1} 'Radio']); % Find button to select
            set(formatPanelHandles.dataTypePanel,'SelectedObject',rightButton); % Ste button group
        else
            set(formatPanelHandles.dataTypePanel,'SelectedObject',[]); % unselected
        end
    end

    function mycloserequest(varargin) %#ok
        % Close figure request turns off data table callback and tidy up.

        if isfield(uitableJavaHandles,'dataTable')
            % Turn off callback to stop it being called
            datatablecallbackcontrol('off');
        end
        delete(fig); % delete figure
        if ~isempty(fid)
            fclose(fid); % close file
        end
    end

    function datatablecallbackcontrol(state)
        % Turn datatable callback on and off to avoid being called lots
        % of times

        if isfield(uitableJavaHandles,'dataTable') && ishandle(uitableJavaHandles.dataTable)
            cm = uitableJavaHandles.dataTable.getTable.getColumnModel; % Get column model object
            switch state
                case 'on'
                    set(cm, 'ColumnSelectionChangedCallback', @dataTableCallback); % Enable
                case 'off'
                    set(cm, 'ColumnSelectionChangedCallback', ''); % Remove/Disable
                otherwise
                    error('unknown state');
            end
        end
    end

    function fileopencallback(varargin) %#ok
        % When 'browse...' button is clicked in preview panel

        cd (pathName); % Change to current or last directory looked at
        currentPathName=pathName;
        currentFileName=fileName;
        [fileName,pathName]= uigetfile({'*.txt;*.csv;*.dat;*.asc;*.ascii';'*.*'},'Text File'); % Get file name

        if isequal(fileName,0) || isequal(pathName,0)
            % Do nothing if fails or exists without chooseing any file
            pathName=currentPathName; % Put back
            fileName=currentFileName; % Put back
        else
            % Pick a file
            set(pushButtonHandles.nextPushButton,'Enable','on');
            numHeaderLines=1;
            columnHeaderRowNumber=1;
            delimiter=',';
            previewpanel;

        end
    end

    function populatepreviewpanel

        % Change title bar
        set(fig,'name','Textscan Tool: Specify Header and Delimiter Info');

        % Set file info uicontrols
        fileInfo=dir([pathName fileName]);
        %fileSizeKB=ceil(fileInfo.bytes/2^10);
        set(previewPanelHandles.fileNameText,'String',[pathName fileName]);

        str=bytes2KMGBytes(fileInfo.bytes);
        set(previewPanelHandles.fileSizeText,'String',['Size: ' str]);

        % Reset headerline info
        set(previewPanelHandles.numHeaderLinesEditBox,'Enable','on');
        set(previewPanelHandles.columnHeadersCheckBox,'Enable','on');
        set(previewPanelHandles.columnHeadersRowEditBox,'Enable','on');
        set(previewPanelHandles.delimiterEditBox,'Enable','on');

        set(previewPanelHandles.numHeaderLinesEditBox,'String',num2str(numHeaderLines));
        set(previewPanelHandles.columnHeadersRowEditBox,'String',num2str(columnHeaderRowNumber));
        set(previewPanelHandles.delimiterEditBox,'String',delimiter);

        % Open file
        fid = fopen([pathName fileName],'r');  % Open text file
        if fid==-1
            error(['File ' fileName ' not found']);
        end

        % Read lines of file as strings
        succeeded=false;
        while ~succeeded
            try
                inputCell=textscan(fid,'%s',numRowsTest,'delimiter','\n','bufsize',bufferSize); % Read strings delimited by a carriage return
                if length(inputCell{1})<100
                    totalNumLinesText=length(inputCell{1});
                end
                succeeded=true;
            catch
                bufferSize=2*bufferSize; % Double buffer size
                fseek(fid,0,'bof'); % Rewind
            end
        end
        lineStrings=inputCell{1};

        % Make strings uitable

        % Adjust column width
        maxLineLength=max(cellfun(@(x) length(x),lineStrings));
        if maxLineLength>200
            columnWidthLineStrings=maxLineLength*5.6;
        end

        uitableJavaHandles.lineStringsTable.setData(lineStrings);
        uitableJavaHandles.lineStringsTable.setColumnWidth(columnWidthLineStrings);
        set(uitableJavaHandles.lineStringsTable,'Visible',1);

    end

    function dataCell2D=textscandatato2Dcell(dataCell)

        % Preallocate dataCell2D
        dataCell2D=cell(length(dataCell{1}),length(dataCell));

        % Convert to a 2D cell
        for k=1:length(dataCell)
            data=dataCell{1,k};
            if iscell(data) % Cell array of strings
                dataCell2D(:,k)=data;
            else % Numeric array
                dataCell2D(:,k)=num2cell(double(data),2); % Fix for unsigned data problem with Java
            end
        end
        dataCell2D(cellfun('isempty',dataCell2D))={'<empty>'};
    end

    function generatefile(varargin) %#ok

        % Write code to file
        str='';
        for line=1:length(codeLines)
            if ~isempty(codeLines{line})
                temp=strrep(codeLines{line},'\','\\'); % Fix escape sequences
                str=[str sprintf([strrep(temp,'%','%%') '\n' ])]; % Fix format
            else
                str=[str sprintf('\n')];
            end
        end

        % Open file for viewing
        com.mathworks.mlservices.MLEditorServices.newDocument(str,true);
        strings = com.mathworks.mlservices.MLEditorServices.builtinGetOpenDocumentNames();
        currName = strings(end);
        com.mathworks.mlservices.MLEditorServices.openDocumentToLine(currName,1);
        
        % If older version
        %if exist('dataimport.m','file')
        %    delete 'dataimport.m'
        %end
        %fid=fopen('dataimport.m','wt');
        %fprintf(fid,'%s',str);
        %fclose(fid);
        %edit('dataimport.m');
    end

    function rowselectcallback(varargin) %#ok

        selectedRow=1+ uitableJavaHandles.lineStringsTable.getTable.getSelectedRow;

        numHeaderLines=selectedRow;
        columnHeaderRowNumber=selectedRow;

        % Get first line of data
        firstDataLineString=lineStrings{selectedRow+1}; % Get first row of data

        % Strip out some characters
        possibleDelimters={' ', ',', ';', '\t'};
        possibleCodes=[32, 44, 59, 9];

        charCodes=double(firstDataLineString'); % Get ASCII codes
        foundDelimiterCodes=...
            (repmat(charCodes,1,length(possibleCodes))==repmat(possibleCodes,length(charCodes),1));
        foundAny=any(foundDelimiterCodes,2);

        foundDelimiters=charCodes(foundAny); % Not numbers

        mostCommonCode=mode(foundDelimiters); % most popular char
        delimiter=possibleDelimters{(mostCommonCode==possibleCodes)'};

        % Set delimiter
        set(previewPanelHandles.numHeaderLinesEditBox,'String',num2str(numHeaderLines));
        set(previewPanelHandles.columnHeadersCheckBox,'Value',1);
        set(previewPanelHandles.columnHeadersRowEditBox,'String',num2str(columnHeaderRowNumber));
        set(previewPanelHandles.delimiterEditBox,'String',delimiter); % find automatically

    end

    function gencodelines(varargin) %#ok
        % Generate M file code lines base on paramters (pop-ups) set

        requestedDataStorageFormatOption=get(importPanelHandles.dataStorageFormatPopup,'Value'); % Get data format option
        numBlockRowsStr=get(importPanelHandles.numRowsEdit,'String'); % Get number of rows from text box
        numBlockRows=str2double(removeCommas(numBlockRowsStr)); % Remove commas is necessary

        startingDataRowStr=get(importPanelHandles.startingDataRowEdit,'String');
        startingDataRow=str2double(removeCommas(startingDataRowStr));
        
        switch requestedDataStorageFormatOption
            case 1 % 1D Cell array of column vectors
                % Generate header and parameters code
                codeLines=[codeSegments('header');...
                    codeSegments('parameters')];
                codeLines=[codeLines;...
                    codeSegments('readBlock')];
                % Put code in list box
                dataStorageFormatOption=requestedDataStorageFormatOption;

                % Calculate size of import variable and time to import
                importToBase=false;
                importDataCallBack % Go ahead and import small piece again in order to calculate size
                importToBase=true;
                dataInfo=whos('data');
                importSize=numBlockRows*(dataInfo.bytes)/numRowsTest; % Scale up
                set(importPanelHandles.importSizeText,'string',['Estimated variable size: ' bytes2KMGBytes(importSize)]);

                importTime=numBlockRows*importTimePiece/numRowsTest; % Estimate import time
                if importTime<1;
                    set(importPanelHandles.importTimeText,'string',['Estimated import time: ' '<1' '(s)']);
                else
                    set(importPanelHandles.importTimeText,'string',['Estimated import time: ' num2str(importTime,'%.2f') '(s)']);
                end

            case 2 % 2D numerical
                % Find unique data types that are not 'ignore'
                formatName=unique(formatNameCell(~strcmp(formatNameCell,'ignore')));
                % Preallocate 2D array
                if length(formatName)==1 && ~any(strcmp(formatNameCell,'string'))...
                        && ~any(strcmp(formatNameCell,'time'));
                    % Generate header and parameters code
                    codeLines=[codeSegments('header');...
                        codeSegments('parameters')];
                    codeLines=[codeLines ;...
                        codeSegments('2DNumerical')];

                    % Calculate size of import variable and time to import
                    importToBase=false;
                    importDataCallBack % Go ahead and import small piece again in order to calculate size
                    importToBase=true;
                    dataInfo=whos('data');
                    importSize=numBlockRows*(dataInfo.bytes)/numRowsTest; % Scale up
                    set(importPanelHandles.importSizeText,'string',['Estimated variable size: ' bytes2KMGBytes(importSize)]);

                    importTime=numBlockRows*importTimePiece/numRowsTest; % Estimate import time
                    if importTime<1;
                        set(importPanelHandles.importTimeText,'string',['Estimated import time: ' '<1' '(s)']);
                    else
                        set(importPanelHandles.importTimeText,'string',['Estimated import time: ' num2str(importTime,'%.2f') '(s)']);
                    end

                else
                    errordlg('Only available for numerical arrays or same type');
                end
                % Put back to old value of 2D Numerical not allowed
                set(importPanelHandles.dataStorageFormatPopup,'Value',dataStorageFormatOption);

            case 3 % 2D Cell
                % Generate header and parameters code
                ignoreColumns=find(strcmp(formatNameCell,'ignore'));
                codeLines=[codeSegments('header');...
                    codeSegments('parameters');...
                    'ignoreColumns=' mat2str(ignoreColumns) ';'];
                codeLines=[codeLines;...
                    codeSegments('2DCell')];
                % Put code in list box
                dataStorageFormatOption=requestedDataStorageFormatOption; % Set it if it reached here

                % Calculate size of import variable and time to import
                importToBase=false;
                importDataCallBack % Go ahead and import small piece again in order to calculate size
                importToBase=true;
                dataInfo=whos('data');
                importSize=numBlockRows*(dataInfo.bytes)/numRowsTest; % Scale up
                set(importPanelHandles.importSizeText,'string',['Estimated variable size: ' bytes2KMGBytes(importSize)]);

                importTime=numBlockRows*importTimePiece/numRowsTest; % Estimate import time
                if importTime<1;
                    set(importPanelHandles.importTimeText,'string',['Estimated import time: ' '<1' '(s)']);
                else
                    set(importPanelHandles.importTimeText,'string',['Estimated import time: ' num2str(importTime,'%.2f') '(s)']);
                end

            otherwise
                errordlg('Unknown data storage option');
        end


    end

    function importDataCallBack(varargin)

        if importToBase
            numBlockRowstoRead=numBlockRows;
            set(importPanelHandles.importDataButton,'String','Busy...');
            drawnow;
        else
            if numBlockRows<numRowsTest
                numBlockRowstoRead=numBlockRows;
            else
                numBlockRowstoRead=numRowsTest;
            end
        end

        dataStorageFormatOption=get(importPanelHandles.dataStorageFormatPopup,'Value');
        fseek(fid,0,'bof'); % Rewind to start
        switch dataStorageFormatOption
            case 1 % 1D Cell array
                % Read data
                tic
                data=textscan(fid,formatString,numBlockRowstoRead,...
                    'headerlines',numHeaderLines+(startingDataRow-1),'delimiter', delimiter);
                importTimePiece=toc;

                % If any date/time columns, then they must be converted to doubles
                found=strcmp(formatStringCell,'%s'); % Find time strings
                if any(found) % Assumes anon function is valid as it has been tested before in GUI
                    stringFnc=get(formatPanelHandles.timeFormatEdit,'String'); % Get anon funtion string for edit box
                    dateNumFncHandle=eval(stringFnc); % Evaluate function Handle from string to create handle
                    data(found)=cellfun(dateNumFncHandle,data(found),'UniformOutput',false); % Convert cell entries
                end

            case 2 % 2D Numerical (all same type)
                if verLessThan('matlab','7.4') % If MATLAB older than verison 7.4
                    tic
                    data=cell2mat(textscan(fid,formatString,numBlockRowstoRead,... % Read data and convert to 2D numerical
                        'headerlines',numHeaderLines+(startingDataRow-1),'delimiter',delimiter));
                    importTimePiece=toc;
                else % More memory efficient using CollectOutput option
                    tic
                    data=textscan(fid,formatString,numBlockRowstoRead,...
                        'headerlines',numHeaderLines+(startingDataRow-1),'delimiter',delimiter,'CollectOutput',true);
                    importTimePiece=toc;
                    data=data{1}; % Extract numerical array
                end

            case 3 % 2D Cell array
                tic
                dataCell=textscan(fid,formatString,numBlockRowstoRead,...  % Read data
                    'headerlines',numHeaderLines+(startingDataRow-1),'delimiter',delimiter);
                importTimePiece=toc;

                % If any date/time columns, then they must be converted to doubles
                found=strcmp(formatStringCell,'%s'); % Find time strings
                if any(found) % Assumes anon function is valid as it has been tested before in GUI
                    stringFnc=get(formatPanelHandles.timeFormatEdit,'String'); % Get anon funtion string for edit box
                    dateNumFncHandle=eval(stringFnc); % Evaluate function Handle from string to create handle
                    dataCell(found)=cellfun(dateNumFncHandle,dataCell(found),'UniformOutput',false); % Convert cell entries
                end

                % Fill up cell
                data=cell(numBlockRowstoRead+1,length(dataCell)); % Preallocate cell
                data(1,:)=columnNames(~strcmp(formatNameCell,'ignore')); % Start with not ignored column names

                for k=1:length(dataCell)
                    if iscell(dataCell{k}) % If cell array of strings...
                        data(2:end,k)=dataCell{k}; %.. just copy in
                    else % If numerical, turn into cells
                        data(2:length(dataCell{k})+1,k)=mat2cell(dataCell{k},ones(length(dataCell{k}),1),1);
                    end
                end

            otherwise
                error('');
        end
        if importToBase
            assignin('base','data',data);
            assignin('base','columnHeaders',columnNames(~strcmp(formatNameCell,'ignore')));         
            clear data; % Clear from GUI workspace
            set(importPanelHandles.importDataButton,'String','Import Data');
        end
    end

    function code=codeSegments(segment)
        switch segment
            case 'header'
                code={...
                    ['function data=' generatedFunction];
                    ['% Import data from ' fileName];
                    ['% Automatically generated ' date];
                    ' ';
                    };
            case 'parameters' % and file open
                numBlocksRowsStr=int2str(numBlockRows);
                code={...
                    '% Define parameters';
                    ['fileName=''' pathName fileName ''';'];
                    ['numHeaderLines=' int2str(numHeaderLines + startingDataRow-1) ';'];                     
                    ['formatString=''' formatString ''';'];
                    ['numRows=' numBlocksRowsStr ';'];
                    ' ';
                    '% Read data from file';
                    'fid=fopen(fileName,''rt'');';
                    };
            case 'readBlock'
                code={... % Start
                    ['data=textscan(fid,formatString,numRows,''headerlines'',numHeaderLines,''delimiter'',''' delimiter ''');'];
                    'status=fclose(fid);';
                    ''};

                found=strcmp(formatStringCell,'%s'); % Find misc strings
                if any(found) % Assumes anon function is valid as it has been tested before in GUI
                    code=[code;... % Optional if any datetime columns
                        '% Convert date/times';
                        ['dateNumFncHandle=' get(formatPanelHandles.timeFormatEdit,'String') ';'];
                        ['data(' mat2str(find(found)) ')=cellfun(dateNumFncHandle,data(' mat2str(find(found)) '),''UniformOutput'',false);'];... %
                        ' ';
                        ];
                end

            case '2DNumerical'
                if verLessThan('matlab','7.4') % If MATLAB older than verison 7.4
                    code=['data=cell2mat(textscan(fid,formatString,numRows,''headerlines'',numHeaderLines,''delimiter'',''' delimiter '''));'];
                else % Use more memory efficient methid using CollectOutput option in 7.4
                    code={
                        ['data=textscan(fid,formatString,numRows,''headerlines'',numHeaderLines,''delimiter'',''' delimiter ''',''CollectOutput'',true); % Option is MATLAB 7.4 feature'];
                        'data=data{1};'}; % Extract numerical array
                end
                code=[code;...
                    'status=fclose(fid);'
                    ' ';
                    ];

            case '2DCell' % With column names
                code={...
                    ['columnHeaderRowNumber=' int2str(columnHeaderRowNumber) ';'];
                    
                    ['columnHeaderRowString=textscan(fid,''%s'',1,''headerlines'',columnHeaderRowNumber-1,''delimiter'',''\n'');'];
                    ['columnHeaderRow=textscan(columnHeaderRowString{1}{1},repmat(''%s'',1,' num2str(numColumns) '),1,''delimiter'',''' delimiter ''');'];
                    ' columnNames=columnHeaderRow;';
                    'columnNames(ignoreColumns)=[];';
                    ['dataCell=textscan(fid,formatString,numRows,''headerlines'',numHeaderLines-columnHeaderRowNumber,''delimiter'',''' delimiter ''');'];
                    ' '};

                % If any time/date string
                found=strcmp(formatStringCell,'%s'); % Find time strings
                if any(found) % Assumes anon function is valid as it has been tested before in GUI
                    code=[code;... % Optional if any datetime columns
                        '% Convert date/times';
                        ['dateNumFncHandle=' get(formatPanelHandles.timeFormatEdit,'String') ';'];
                        ['dataCell(' mat2str(find(found)) ')=cellfun(dateNumFncHandle,dataCell(' mat2str(find(found)) '),''UniformOutput'',false);'];... %
                        ' ';
                        ];
                end

                code=[code;... % Ending
                    {...
                    '% Convert to 2D Cell'
                    'data=cell(numRows+1,length(dataCell)); %Preallocate';
                    'data(1,:)=columnNames; % Start with column names';
                    '';
                    'for k=1:length(dataCell)';
                    '    if iscell(dataCell{k}) % If cell array of strings';
                    '        data(2:end,k)=dataCell{k};';
                    '    else % If numerical';
                    '       data(2:end,k)=mat2cell(dataCell{k},ones(numRows,1),1);';
                    '   end';
                    'end';
                    'status=fclose(fid);';
                    };
                    ];

            otherwise
                error('Unknown code segment');
        end
    end
    function countLinesCallBack(varargin)
        % Count the number of lines in the text file

        if ~isempty(fileName) % As long as there is a file opened
            if isempty(totalNumLinesText) % If hasn't been calculate yet
                fileInfo=dir([pathName fileName]); % File stats
                fileSize=ceil(fileInfo.bytes); % File size
                h = waitbar(0,'Counting Lines'); % Create wait bar
                set(h,'name','Counting Lines in Text File'); % Set figure name
                set(h,'Color',get(0,'DefaultUipanelBackgroundColor'));
                amountToRead=1e5; % Amount of bytes to read in at a time in one chunk
                numBlocks=ceil(fileSize/amountToRead); % NUmber of chunks to read

                %% Read
                totalNumLinesText=1; % Initialize count of number of lines
                fseek(fid,0,'bof'); % Rewind to start
                for k=1:numBlocks % ~feof(fid)
                    l=fread(fid,amountToRead,'*uint8'); % Read chunk of uint8s
                    totalNumLinesText=totalNumLinesText+length(find(l==10)); % Find line feeds and count
                    percent=(k*amountToRead)/fileSize; % Percent read in
                    waitbar(percent,h,[int2str(100*percent) '% Complete']); % Update wait bar
                end
                close(h); % Delete waitbar
            end
            numBlockRows=totalNumLinesText-numHeaderLines;
            set(importPanelHandles.totalNumLinesText,'String',['Total No. lines: ' numbercommaformat(int2str(totalNumLinesText))]); % Total lines
            set(importPanelHandles.numRowsEdit,'String',numbercommaformat(int2str(totalNumLinesText-numHeaderLines))); % Less Header lines
            gencodelines
        end
    end

    function output=numbercommaformat(input) %#ok
        % Converts a number string to one with commas seperating thousands
        % Example
        %   numbercommaformat('1088880')
        %   ans =
        %   1,088,880

        input=fliplr(input); % Reverse string
        count=1; % Start of counter to output string
        for k=1:length(input) % For each character in inpout string
            output(count)=input(k); %#ok Copy over to new string
            count=count+1; % Increment counter
            if mod(k,3)==0 % if devisible by 3
                output(count)=','; %#ok Add a comma
                count=count+1; % Increment counter
            end
        end
        if strcmp(output(end),',') % If comma at the end
            output(end)=[]; % remove it
        end
        output=fliplr(output); % Reverse string back

    end

    function str=bytes2KMGBytes(bytes)
        % Returns string for number in bytes, KB, MB of GB
        if bytes<2^10 % Bytes
            str=[ num2str(bytes) ' B'];
        elseif bytes>=2^10 && bytes<2^20
            str=[ num2str(bytes/2^10,'%.2f') ' KB'];
        elseif bytes>=2^20 && bytes<2^30
            str=[ num2str(bytes/2^20,'%.2f') ' MB'];
        else
            str=[ num2str(bytes/2^30,'%.2f') ' GB'];
        end
    end

    function numberStr=removeCommas(numberStrWithCommas)
        % Remove any commas from a big number string
        numberStr=numberStrWithCommas(numberStrWithCommas~=',');
    end

    function helpcallback(varargin)
        % Help button callback
        web(['file:' which('textscantoolhelp.html')])
    end

end
