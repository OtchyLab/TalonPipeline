
function varargout = Metermeter(varargin)
% METERMETER MATLAB code for Metermeter.fig
%      Last modified by TMO June 28, 2014
%
%      METERMETER, by itself, creates a new METERMETER or raises the existing
%      singleton*.
%
%      H = METERMETER returns the handle to a new METERMETER or the handle to
%      the existing singleton*.
%
%      METERMETER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in METERMETER.M with the given input arguments.
%
%      METERMETER('Property','Value',...) creates a new METERMETER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Metermeter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Metermeter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Metermeter

% Last Modified by GUIDE v2.5 05-Jun-2015 13:53:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Metermeter_OpeningFcn, ...
                   'gui_OutputFcn',  @Metermeter_OutputFcn, ...
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


% --- Executes just before Metermeter is made visible.
function Metermeter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Metermeter (see VARARGIN)

%Color definitions for iterative plotting
handles.colorcell = {'b', 'r','g', 'k', 'm', 'c', 'y', 'b', 'r','g', 'k', 'm', 'c', 'y'};
handles.symbolcell ={'s','o', 'x',  '+', '*', 's','o', 'x',  '+', '*', 's','o', 'x',  '+', '*'};



% Choose default command line output for Metermeter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Metermeter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Common functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dateVect,hourVect] = extractTimeVects(filenames)
%Extract the date and hour from the filenames
if iscellstr(filenames)
    rendNum = size(filenames,2);
else
    rendNum = size(filenames,1);
end

for i = 1:rendNum
    sp{i} = regexp(filenames{i},'_','split');
    dateVect(i) = datenum([sp{i}{3} '-' sp{i}{4} '-' sp{i}{5}], 'yyyy-mm-dd');
    hourVect(i) = str2num(sp{i}{6}) + ((str2num(sp{i}{7}) + (str2num(sp{i}{8}(1:2))/60)) /60);
end

function IntMat = extractIntervals(p, q, templatesyllBreaks)
%Geat template break points
a = transpose(templatesyllBreaks);
a = a(:);

%Interval Durations from the paths and template files
rendNum = length(p);
IntMat = [];
for i = 1:rendNum
    %Rover rendition path and rendition length
    path = [p{i},q{i}];
    
    %Calculate intervals and add to the stack
    [warpedOut] = getWarpedStarts(path,a(:));
    IntMat = [IntMat; diff(warpedOut)'];
end

function IntMat = calcIntervals(p, q, templatesyllBreaks)
%Extract all interval durations from the paths and templatebreaks

%Get the template syllable breaks from the passed file
a = transpose(templatesyllBreaks);
a = a(:);

%Interval Durations from the paths and template files
rendNum = length(p);
IntMat = [];
path = {};
warpedOut = {};

parfor i = 1:rendNum
    %Rover rendition path and rendition length
    path{i} = [p{i},q{i}];
    
    %Calculate intervals and add to the stack
    warpedOut{i} = getWarpedStarts(path{i},a(:));
    IntMat(i,:) = diff(warpedOut{i})';
end

function subIndices = parseIndices(handles, idx)
%Get parsing paramters from GUI controls
parseRange = str2num(get(handles.edit_parseTimes,'String'));
parseNum = str2num(get(handles.edit_numMotifs,'String'));

%Cycle through each bin and generate file index
[numBins, ~] = size(parseRange);
subIndices = [];
for i = 1:numBins;
    %Raw index
    Indx = find(handles.dataset(idx).timeVect >= parseRange(i,1) & handles.dataset(idx).timeVect <= parseRange(i,2));
    
    if ~isempty(Indx)
        %Reduce to desired rend num
        if parseNum == -1
            pIndx = Indx;
        else
            if length(Indx) > parseNum
                pIndx = Indx(1:parseNum);
            else
                pIndx = Indx;
            end
        end
        pnts = [pIndx(1), pIndx(end)];
    else
        pnts = [];
    end
    
    %Add to the stack
    subIndices = [subIndices; pnts];
end

function pData = parseData(handles)

%Determine the number of datasets in list
if iscell(handles.datasetFilename)
    numDatasets = size(handles.datasetFilename,2);
else
    numDatasets = size(handles.datasetFilename,1);
end

%Parse out the datasets
s = 1;
pData = [];
for i = 1:numDatasets
    %Generate parsing indices according to user selections
    subIndices = parseIndices(handles, i);

    %Unpack the parsed data using the 
    for j = 1:size(subIndices,1);
        %Generate the index
        idx = subIndices(j,1):subIndices(j,2);
        
        %Apply index
        pData(s).intervals = handles.dataset(i).intervals(idx,:);
        pData(s).filenames = handles.dataset(i).renditionFilenames(idx);
        pData(s).date = handles.dataset(i).dateVect(idx);
        pData(s).time = handles.dataset(i).timeVect(idx);
        
        %Update the pointer
        s = s + 1;
    end
end

function parseNames = updateParseNames(handles)
%Determine the number of datasets in list
if iscell(handles.datasetFilename)
    numDatasets = size(handles.datasetFilename,2);
else
    numDatasets = size(handles.datasetFilename,1);
end

%Parse out the datasets
s = 1;
parseNames = [];
tm = str2num(get(handles.edit_parseTimes,'String'));
for i = 1:numDatasets
    sp = regexp(handles.datatitles{i}, '_', 'split');
    nm = [sp{1} '_' sp{2}];
    
    %Unpack the parsed data using the 
    for j = 1:size(tm,1);
        %Set Name       
        parseNames{s} = [nm ' @ ' num2str(tm(j,1))];
        
        %Update the pointer
        s = s + 1;
    end
end

function dataOut = normData(ind, dataIn)
[m, n] = size(dataIn);

for i = 1:n
    dataOut(:,i) = dataIn(:,i)./mean(dataIn(ind,i));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parsing functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function list_SylGap_Callback(hObject, eventdata, handles)

function list_SylGap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function radio_SylsGaps_Callback(hObject, eventdata, handles)
%Checking this box, sets the paired box to false
set(handles.radio_SylNGap,'Value',0);

%You may not uncheck the box
if ~get(handles.radio_SylsGaps,'Value')
    set(handles.radio_SylsGaps,'Value',1);
end

%Update the listbox values
if isfield(handles, 'sgLabels')
    set(handles.list_SylGap, 'String', [handles.sgLabels; 'Total'], 'Value', 1)
end

guidata(hObject, handles);

function radio_SylNGap_Callback(hObject, eventdata, handles)
%Checking this box, sets the paired box to false
set(handles.radio_SylsGaps,'Value',0);

%You may not uncheck the box
if ~get(handles.radio_SylNGap,'Value')
    set(handles.radio_SylNGap,'Value',1);
end

%Update the listbox values
if isfield(handles, 'sgLabels')
    set(handles.list_SylGap, 'String', [handles.sNgLabels; 'Total'], 'Value', 1)
end

guidata(hObject, handles);

function edit_parseTimes_Callback(hObject, eventdata, handles)
%Update parse names
handles.parseNames = updateParseNames(handles);
set(handles.list_norm, 'String', handles.parseNames);

guidata(hObject, handles);

function edit_parseTimes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_numMotifs_Callback(hObject, eventdata, handles)

function edit_numMotifs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Data load functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_addDataset_Callback(hObject, eventdata, handles)
%This function adds new dataset files to the filelist, sorts them, and
%displays in the listbox. Some basic processing is done, but no data is
%loaded from the source. (That requires an additional button click.)

%Get the locations of the processed files (multiselect is possible)
[temp,path] = uigetfile('*.mat','Select the dataset file.','MultiSelect','on');

if isequal(temp,0) || isequal(path,0)
    %If they hit cancel or choose something dumb, generate error and end.
    set(handles.text_message,'String','Dataset location invalid or cancelled. Pick a valid file.');
else
    %If they don't, then parse the locations.
    fnames = [];
    if iscell(temp)
        %If there is more than one files selected, struct it
         for i = 1:size(temp,2)
             fnames{i} = char(temp{i});
             pathTot{i} = [path char(temp{i})];
         end
    else
        %Otherwise, just copy it out.
        fnames = temp;
        pathTot = [path,temp];
        i = 1;
    end

    if ~isfield(handles, 'datasetFilename')
        handles.datasetFilename = [];
    end
    
     if ~isempty(handles.datasetFilename)
%          button = questdlg('Do you want to add this dataset(s) to the current list or replace them?','Add an annotation','Add to current','Replace current','Add to current');
%          if (strcmp(button,'Replace current'))
%             %Add here whatever other shit needs to be cleared from memory
%             if isfield(handles,'dataset')
%                 handles = rmfield(handles,'dataset');
%             end
%             if isfield(handles,'procSSData')
%                 handles = rmfield(handles,'procSSData');
%             end
%             if isfield(handles, 'procCount')
%                 handles = rmfield(handles, 'procCount');
%             end
%             if isfield(handles,'long')
%                 handles = rmfield(handles,'long');
%             end
%             
%             %Copy over the new files
%              handles.datasetFilename = pathTot;
%              handles.datatitles = fnames;
% 
%              %Sort the datatitles to alphabetically (which is also temporally if 
%              %only one bird) and use the index to sort the datasetFilename as well
%              if iscell(handles.datatitles)
%                 [handles.datatitles, ind] = sort(handles.datatitles);
%                 handles.datasetFilename = handles.datasetFilename(ind);
%              end
% 
%              %Set the message to show annotation replaced
%              set(handles.listbox_datasets,'String',handles.datatitles)
%              set(handles.text_message,'String',['Added ' num2str(i) ' new datasets to the active list']);
%              
%          %This is the code for adding to the existing list 
%          elseif (strcmp(button,'Add to current'))
%              %Check if there are multiple files in queue already
%              if iscell(handles.datasetFilename)
%                  dFtemp = handles.datasetFilename{:};
%                  dTtemp = handles.datatitles{:};
%              else
%                  dFtemp = handles.datasetFilename;
%                  dTtemp = handles.datatitles;
%              end
%              
%              %Check if there are multiple files to add
%              if iscell(fnames)
%                 handles.datasetFilename = {dFtemp pathTot{:}}';
%                 handles.datatitles = {dTtemp fnames{:}}';
%              else
%                 handles.datasetFilename = {dFtemp pathTot}';
%                 handles.datatitles = {dTtemp fnames}';                 
%              end
%              
%              %Remove duplicates
%              handles.datatitles = unique(handles.datatitles);
%              handles.datasetFilename = unique(handles.datasetFilename);
%                          
%              %Sort the datatitles to alphabetically (which is also temporally if 
%              %only one bird) and use the index to sort the datasetFilename as well
%              if iscell(handles.datatitles)
%                 [handles.datatitles, ind] = sort(handles.datatitles);
%                 handles.datasetFilename = handles.datasetFilename(ind);
%              end
% 
%              %Set the message to show annotation replaced
%              set(handles.listbox_datasets,'String',handles.datatitles)
%              set(handles.text_message,'String',['Added ' num2str(i) ' new datasets to the active list']);
%          end
     %This is the code for loading the first/only folder    
     else
         handles.datasetFilename = pathTot;
         datatitles = fnames;

         %Sort the datatitles alphabetically (which is also chronologically if 
         %only one bird) and use the index to sort the datasetFilename as well
         if iscell(datatitles)
            [handles.datatitles, ind] = sort(datatitles);
            handles.datasetFilename = handles.datasetFilename(ind);
         else
            handles.datatitles = datatitles;
            handles.datasetFilename = handles.datasetFilename;
         end

         %Set the message to show annotation replaced
         set(handles.listbox_datasets,'String',handles.datatitles)
         set(handles.text_message,'String',['Added ' num2str(i) ' new datasets to the active list']);
     end
     
end
guidata(hObject, handles);

function push_clearDatasets_Callback(hObject, eventdata, handles)
%Get user confirmation on clearing the list
button = questdlg('Are you sure you want to clear all datasets?','Clear datasets?','Clear''em All!','Nooooo!','Nooooo!');

if (strcmp(button,'Clear''em All!'))
    %Clear pointers to the datasets and their labels
     handles.datasetFilename = [];
     handles.datatitles = [];
     
    %Add here whatever other shit needs to be cleared from memory


    %Set the message to show datasets replaced
    set(handles.listbox_datasets,'String',handles.datatitles)
    set(handles.text_message,'String','All datasets cleared.');
else
    set(handles.text_message,'String','Datasets are unchanged.');
end

guidata(hObject, handles);

function push_loadDataset_Callback(hObject, eventdata, handles)
%Load the data from the saved .dat and .wav files to the workspace
if isempty(handles.datasetFilename)
    set(handles.text_message,'String','No datasets found to load. Check list, add files and try again.')
else
     %Work through each file to extract the interval data
    if iscell(handles.datasetFilename)
        numDatasets = size(handles.datasetFilename,2);
    else
        numDatasets = size(handles.datasetFilename,1);
    end

    %Cycle through each dataset, load from file, and extract needed info.
    h = waitbar(0,'Loading datasets and extracting intervals. Please wait...');
    handles.dataset = [];
    for i = 1:numDatasets
        %Load specific variables from file
        if iscell(handles.datasetFilename)
            load(handles.datasetFilename{i},'data','filenames','sequence')
        else
            load(handles.datasetFilename,'data','filenames','sequence')
        end
        
        %Parse the filenames for the time of recordings
         [dateVect,timeVect] = extractTimeVects(filenames);
        
        %Parse the sequence data to strip out the rendition count
        sp = regexp(sequence,' ','split');
        seqCheck(i) = str2double(sp(1));

        % Generate the interval labels for the first dataset loaded... this could cause problems if datasets aren't for the
        % same sequence, but this would ikely show up somewhewre else and vcrash the whole thing...
        if i ==1;
            seq = str2double(sp(1));
            %Create Interval Labels
            sgLabels = {};
            sNgLabels = {};
            sLabels = {};
            for j = 1:length(sp{1})
                %Labels for breaking out gaps from syllables
                sgLabels = [sgLabels; sp{1}(j)];
                if j ~= length(sp{1}) 
                    sgLabels = [sgLabels; [sp{1}(j) 'G']];
                end

                %Labels for joinging syllables and gaps
                if j ~= length(sp{1})
                    sNgLabels = [sNgLabels; [sp{1}(j) '+G']];
                else
                    sNgLabels = [sNgLabels; sp{1}(j)];
                end

                %Labels for just syllables
                sLabels = [sLabels; sp{1}(j)];
            end
            %Copy out
            handles.sgLabels = sgLabels;
            handles.sNgLabels = sNgLabels;
            handles.sLabels = sLabels;
        end
        
        %Extract S/G intervals for all renditions
        IntMat = extractIntervals(data.p, data.q, data.templatesyllBreaks);
        
        %Copy out data to the handles structure
        handles.dataset(i).renditionFilenames = filenames;
        handles.dataset(i).templatesyllBreaks = data.templatesyllBreaks;
        handles.dataset(i).p = data.p;
        handles.dataset(i).q = data.q;
        handles.dataset(i).intervals = IntMat;
        handles.dataset(i).dateVect = dateVect;
        handles.dataset(i).timeVect = timeVect;

        %Clear out already processed data to free memory
        clear('data','filenames','sequence');
        waitbar(i/numDatasets)
    end
    set(handles.text_message,'String','Dataset creation now completed.')
    close(h) 
end

%Update parse names
handles.parseNames = updateParseNames(handles);
set(handles.list_norm, 'String', handles.parseNames);

%Perhaps should check to make sure all of the sequences are the same
if length(unique(seqCheck)) > 1
    warndlg('Not all datasets appear to have the same syllable sequence. Likely problems ahead.');
    uiwait;
end

guidata(hObject, handles);

function listbox_datasets_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_datasets_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Save/export buttons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_plotSimpleMean_Callback(hObject, eventdata, handles)
%Simple mean +/- std plot

%Parse Data
handles.pData = parseData(handles);

%Get interval selections from GUI controls
vals = get(handles.list_SylGap, 'Value');
contents = get(handles.list_SylGap, 'String');
if get(handles.radio_SylsGaps,'Value')
    intType = 1;
else
    intType = 0;
end

%Do the stats on the parsed data
numBins = length(handles.pData);
for i = 1:numBins
    for j = 1:length(vals)
        if vals(j) ~= length(contents) && intType %Syl and Gap, not total, selected
            int_m(i, j) = mean(handles.pData(i).intervals(:,vals(j)));
            int_std(i, j) = std(handles.pData(i).intervals(:,vals(j)));
        elseif vals(j) ~= length(contents) && ~intType %Syl+Gap, not total, selected
            if vals(j) ~= (length(contents) -1) %requires summing syl and gap
                pntr = (2*vals(j)-1):(2*vals(j));
            else
                pntr = (2*vals(j)-1); %last syl, so no gap to add
            end
            int_m(i, j) = mean(sum(handles.pData(i).intervals(:,pntr), 2));
            int_std(i, j) = std(sum(handles.pData(i).intervals(:,pntr), 2));
        elseif vals(j) == length(contents) %Total selected
            int_m(i, j) = mean(sum(handles.pData(i).intervals, 2));
            int_std(i, j) = std(sum(handles.pData(i).intervals, 2));
        end
    end
    
    %Calc/format timepoints
    dt(i) = (handles.pData(i).date(1) + (handles.pData(i).time(1)/24));
end

%Normalize data, if selected
if get(handles.check_norm, 'Value')
    int_std = int_std./int_m;
    int_m = normData(get(handles.list_norm,'Value'), int_m);
end

%Get zeroDay for plotting
zeroDay = datenum(get(handles.edit_zeroDay, 'String'), 'mm/dd/yy');

%Plot output
if get(handles.check_exportPlot, 'Value')
    figure;
else
    axes(handles.axes1);
end
cla; hold on
for i = 1:length(vals)
    errorbar(dt-zeroDay, int_m(:,i), int_std(:,i), [':' handles.symbolcell{i} handles.colorcell{i}], 'LineWidth', 1.5)
end
axis auto
set(gca, 'Box', 'off', 'TickDir', 'out')
xlabel('Time (Days)', 'FontSize', 20)
if get(handles.check_norm, 'Value')
    ylabel('Norm Duration (%)', 'FontSize', 20)
else
    ylabel('Duration (ms)', 'FontSize', 20)
end
set(gca, 'LineWidth', 3, 'FontSize', 20)
title('Mean Motif Length', 'FontSize', 20);
legend(contents(vals)); legend('boxoff')
hold off

%Copy useful data out to save structure
handles.saveData = [];
handles.saveData.time = dt-zeroDay;
handles.saveData.int_m = int_m;
handles.saveData.int_std = int_std;

handles.saveData.dataType = contents(vals);
handles.saveData.parseTime = str2num(get(handles.edit_parseTimes,'String'));
handles.saveData.parseNum = str2num(get(handles.edit_numMotifs,'String'));
handles.saveData.normed = get(handles.check_norm, 'Value');
handles.saveData.normRef = get(handles.list_norm,'Value');
handles.saveData.fileSets = get(handles.list_norm,'String');

guidata(hObject, handles);

function push_saveRaw_Callback(hObject, eventdata, handles)

function push_saveParsed_Callback(hObject, eventdata, handles)
dataParsed = handles.pData;

[SaveName,SavePath] = uiputfile('*.mat');
save([SavePath SaveName],'dataParsed');

function push_saveFinal_Callback(hObject, eventdata, handles)
dataOut = handles.saveData;

[SaveName,SavePath] = uiputfile('*.mat');
save([SavePath SaveName],'dataOut');



function push_plotAll_Callback(hObject, eventdata, handles)
%Plot all of the intervals as scatter plot

%Parse Data
handles.pData = parseData(handles);

%Get interval selections from GUI controls
vals = get(handles.list_SylGap, 'Value');
contents = get(handles.list_SylGap, 'String');
if get(handles.radio_SylsGaps,'Value')
    intType = 1;
else
    intType = 0;
end

%Do the stats on the parsed data
numBins = length(handles.pData);
int = [];
dt = [];
for i = 1:numBins
    len = size(handles.pData(i).intervals,1);
    indx{i} = (size(int,1)+1):(size(int,1)+len);
    for j = 1:length(vals)
        if vals(j) ~= length(contents) && intType %Syl and Gap, not total, selected
            int(indx{i}, j) = handles.pData(i).intervals(:,vals(j));
        elseif vals(j) ~= length(contents) && ~intType %Syl+Gap, not total, selected
            if vals(j) ~= (length(contents) -1) %requires summing syl and gap
                pntr = (2*vals(j)-1):(2*vals(j));
            else
                pntr = (2*vals(j)-1); %last syl, so no gap to add
            end
            int(indx{i}, j) = sum(handles.pData(i).intervals(:,pntr), 2);
        elseif vals(j) == length(contents) %Total selected
            int(indx{i}, j) = sum(handles.pData(i).intervals, 2);
        end
    end
    
    %Calc/format timepoints
    dt(indx{i}) = (handles.pData(i).date + (handles.pData(i).time/24));
end

%Normalize data, if selected
if get(handles.check_norm, 'Value')
    mInd = [];
    for i = get(handles.list_norm,'Value')
        mInd = [mInd indx{i}];
    end
    int = normData(mInd, int);
end

%Get zeroDay for plotting
zeroDay = datenum(get(handles.edit_zeroDay, 'String'), 'mm/dd/yy');

%Plot output
if get(handles.check_exportPlot, 'Value')
    figure;
else
    axes(handles.axes1);
end
cla; hold on
for i = 1:length(vals)
    scatter(dt-zeroDay, int(:,i), ['.' handles.colorcell{i}])
end
axis auto
set(gca, 'Box', 'off', 'TickDir', 'out')
xlabel('Time (Days)', 'FontSize', 20)
if get(handles.check_norm, 'Value')
    ylabel('Norm Duration (%)', 'FontSize', 20)
else
    ylabel('Duration (ms)', 'FontSize', 20)
end
set(gca, 'LineWidth', 3, 'FontSize', 20)
title('Motif Length', 'FontSize', 20);
legend(contents(vals)); legend('boxoff')
hold off

%Copy useful data out to save structure
handles.saveData = [];
handles.saveData.time = dt-zeroDay;
handles.saveData.int = int;

handles.saveData.dataType = contents(vals);
handles.saveData.parseTime = str2num(get(handles.edit_parseTimes,'String'));
handles.saveData.parseNum = str2num(get(handles.edit_numMotifs,'String'));
handles.saveData.normed = get(handles.check_norm, 'Value');
handles.saveData.normRef = get(handles.list_norm,'Value');
handles.saveData.fileSets = get(handles.list_norm,'String');

guidata(hObject, handles);

function push_plotSimpCV_Callback(hObject, eventdata, handles)
%Simple CV

%Parse Data
handles.pData = parseData(handles);

%Get interval selections from GUI controls
vals = get(handles.list_SylGap, 'Value');
contents = get(handles.list_SylGap, 'String');
if get(handles.radio_SylsGaps,'Value')
    intType = 1;
else
    intType = 0;
end

%Do the stats on the parsed data
numBins = length(handles.pData);
for i = 1:numBins
    for j = 1:length(vals)
        if vals(j) ~= length(contents) && intType %Syl and Gap, not total, selected
            int_m(i, j) = mean(handles.pData(i).intervals(:,vals(j)));
            int_std(i, j) = std(handles.pData(i).intervals(:,vals(j)));
        elseif vals(j) ~= length(contents) && ~intType %Syl+Gap, not total, selected
            if vals(j) ~= (length(contents) -1) %requires summing syl and gap
                pntr = (2*vals(j)-1):(2*vals(j));
            else
                pntr = (2*vals(j)-1); %last syl, so no gap to add
            end
            int_m(i, j) = mean(sum(handles.pData(i).intervals(:,pntr), 2));
            int_std(i, j) = std(sum(handles.pData(i).intervals(:,pntr), 2));
        elseif vals(j) == length(contents) %Total selected
            int_m(i, j) = mean(sum(handles.pData(i).intervals, 2));
            int_std(i, j) = std(sum(handles.pData(i).intervals, 2));
        end
    end
    
    %Calc/format timepoints
    dt(i) = (handles.pData(i).date(1) + (handles.pData(i).time(1)/24));
end

%Calc the CV
int_cv = int_std./int_m;

%Normalize data, if selected
if get(handles.check_norm, 'Value')
    int_cv = normData(get(handles.list_norm,'Value'), int_cv);
end

%Get zeroDay for plotting
zeroDay = datenum(get(handles.edit_zeroDay, 'String'), 'mm/dd/yy');

%Plot output
if get(handles.check_exportPlot, 'Value')
    figure;
else
    axes(handles.axes1);
end
cla; hold on
for i = 1:length(vals)
    plot(dt-zeroDay, int_cv(:,i), [':' handles.symbolcell{i} handles.colorcell{i}])
end
axis auto
set(gca, 'Box', 'off', 'TickDir', 'out')
xlabel('Day & Time')
if get(handles.check_norm, 'Value')
    ylabel('Norm CV (%)')
else
    ylabel('CV (ms)')
end
legend(contents(vals)); legend('boxoff')
hold off

%Copy useful data out to save structure
handles.saveData = [];
handles.saveData.time = dt-zeroDay;
handles.saveData.int_cv = int_cv;

handles.saveData.dataType = contents(vals);
handles.saveData.parseTime = str2num(get(handles.edit_parseTimes,'String'));
handles.saveData.parseNum = str2num(get(handles.edit_numMotifs,'String'));
handles.saveData.normed = get(handles.check_norm, 'Value');
handles.saveData.normRef = get(handles.list_norm,'Value');
handles.saveData.fileSets = get(handles.list_norm,'String');

guidata(hObject, handles);

function edit_zeroDay_Callback(hObject, eventdata, handles)
a=1

function edit_zeroDay_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function list_norm_Callback(hObject, eventdata, handles)

function list_norm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function check_norm_Callback(hObject, eventdata, handles)

function check_exportPlot_Callback(hObject, eventdata, handles)

function push_varDecomp_Callback(hObject, eventdata, handles)

function push_button1_Callback(hObject, eventdata, handles)


function check_sumInts_Callback(hObject, eventdata, handles)

