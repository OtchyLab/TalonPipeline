function varargout = Talon(varargin)
% TALON MATLAB code for Talon.fig
%      TALON, by itself, creates a new TALON or raises the existing
%      singleton*.
%
%      H = TALON returns the handle to a new TALON or the handle to
%      the existing singleton*.
%
%      TALON('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TALON.M with the given input arguments.
%
%      TALON('Property','Value',...) creates a new TALON or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Talon_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Talon_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Talon

% Last Modified by GUIDE v2.5 08-Feb-2016 18:00:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Talon_OpeningFcn, ...
                   'gui_OutputFcn',  @Talon_OutputFcn, ...
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

function Talon_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

%Set up axes 
set(handles.axes_raster, 'ButtonDownFcn', @cb_audio_click,'XTick',[],'YTick',[]);
set(handles.axes_stats, 'XTick', [], 'YTick', []);
set(handles.axes_template,'XTick',[],'YTick',[]);
set(handles.axes_aux, 'XTick', [], 'YTick', []);

set(handles.popup_alignType,'Value',1);

handles.fs=44150;
handles.datasetFilename =[];
handles.datatitles = [];
handles.celltitles = [];
handles.rastered = false;

%Plot colors
t = get(0,'defaultAxesColorOrder');
handles.colors = [t; rand(23,3)]; %30 initial plotting colors to work with

% Update handles structure
guidata(hObject, handles);

function varargout = Talon_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Independent, in-GUI functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cb_click(hObject, evnt)
handles = guidata(hObject);

%Get the current click location, relative to the axes
click = get(gca,'CurrentPoint');

%If double click then show the record number at click  
handles.spotNum = ceil(click(1,2));

%Update marker and spotlight
handles = updateSpotlight(handles);

display(['Got it: ' num2str(handles.spotNum)])
    
guidata(hObject, handles);

function syncData = crossRef(handles)
%This function matches motifs to cell file entries using the filenum as the matching identifer. It maintains the grouping of data by cells

%Initialize variables
syncData = [];
preLag = str2double(get(handles.edit_lag, 'String'))/1000; %premotor song lag for parsing spike times

%Cycle through each cell
numCells = length(handles.cell);
for i = 1:numCells
    %Retrieve filenum field vectors
    cFilenums = getFieldVector(handles.cell{i}, 'filenum');
    
    %Reset pointer value and temp output
    pntr = 1;
    curSet = [];
    %Cycle through each of the cell filenums
    numCellFiles = length(cFilenums);
    for j = 1:numCellFiles
        %For each, locate the matching indices in the dataset (can be more than one)
        motifIdx = find(handles.filenums == cFilenums(j));
        
        %Check to make sure the match was found
        if ~isempty(motifIdx)
            %Create composite record for each of the found indices
            for k = 1:length(motifIdx)
                %Fields defined in the dataset
                curSet(pntr).filename = handles.filenames{motifIdx(k)};
                curSet(pntr).filenum = handles.filenums(motifIdx(k));
                curSet(pntr).datenum = getFileTime(handles.filenames{motifIdx(k)})+(handles.snipTimes(motifIdx(k),1)/(60*60*24)); %Recording time, plus the offset to motif start
                curSet(pntr).snipTimes = handles.snipTimes(motifIdx(k),:);
                curSet(pntr).pFull = handles.pFull{motifIdx(k)};
                curSet(pntr).drugStatus = handles.drugStatus{motifIdx(k)};
                curSet(pntr).directStatus = handles.directStatus{motifIdx(k)};
                
                %Fields defined in the cell file
                curSet(pntr).cellname = handles.celltitles{i};
                curSet(pntr).channel = handles.cell{i}(j).channel;
                
                %The spike selections account for both the snipTime offset (from the original file) and the pre-motor lead (specified in the GUI; ~15ms for RA)
                sp = handles.cell{i}(j).spikes; noi = handles.cell{i}(j).noise; offset = handles.snipTimes(motifIdx(k),1); %to improve readability below
                curSet(pntr).spikes = sp(sp>=(handles.snipTimes(motifIdx(k),1)-preLag) & sp<=(handles.snipTimes(motifIdx(k),2)-preLag)) - offset; 
                curSet(pntr).noise = noi(noi>=(handles.snipTimes(motifIdx(k),1)-preLag) & noi<=(handles.snipTimes(motifIdx(k),2)-preLag)) - offset;
                
                %Increment pointer value
                pntr = pntr + 1;
            end
        else
            %If there is no matching number, throw an error/message and continue on
            display(['Yo! I''m missing a matching motif record for cell file ' num2str(cFilenums(j)) '. Check it out!'])
        end
        
    end
    
    %Push the local variable to the output variable
    syncData{i} = curSet;
end

function syncData = crossRef2(handles)
%This function matches motifs to cell file entries using the filenum as the matching identifer.
%All matches are grouped in the same structure.

%Initialize variables
syncData = [];
curSet = [];
preLag = str2double(get(handles.edit_lag, 'String'))/1000; %premotor song lag for parsing spike times
pntr = 1;

%Cycle through each cell
numCells = length(handles.cell);
for i = 1:numCells
    %Retrieve filenum field vectors
    cFilenums = getFieldVector(handles.cell{i}, 'filenum');
    
    %Cycle through each of the cell filenums
    numCellFiles = length(cFilenums);
    for j = 1:numCellFiles
        %For each, locate the matching indices in the dataset (can be more than one)
        motifIdx = find(handles.filenums == cFilenums(j));
        
        %Check to make sure the match was found
        if ~isempty(motifIdx)
            %Create composite record for each of the found indices
            for k = 1:length(motifIdx)
                %Sync this file if prereqs are met
                if ~isempty(handles.cell{i}(j).spikes) && ~handles.dataFlag(motifIdx(k))
                    %Fields defined in the dataset
                    curSet(pntr).filename = handles.filenames{motifIdx(k)};
                    curSet(pntr).filenum = handles.filenums(motifIdx(k));
                    curSet(pntr).datenum = getFileTime(handles.filenames{motifIdx(k)})+(handles.snipTimes(motifIdx(k),1)/(60*60*24)); %Recording time, plus the offset to motif start
                    curSet(pntr).snipTimes = handles.snipTimes(motifIdx(k),:);
                    curSet(pntr).pFull = handles.pFull{motifIdx(k)};
                    curSet(pntr).drugStatus = handles.drugStatus{motifIdx(k)};
                    curSet(pntr).directStatus = handles.directStatus{motifIdx(k)};
                    
                    %Fields defined in the cell file
                    curSet(pntr).cellname = handles.celltitles{i}(1:end-4);
                    curSet(pntr).channel = handles.cell{i}(j).channel;
                    
                    %The spike selections account for both the snipTime offset (from the original file) and the pre-motor lead (specified in the GUI; ~15ms for RA)
                    sp = handles.cell{i}(j).spikes;
                    noi = handles.cell{i}(j).noise;
                    offset = handles.snipTimes(motifIdx(k),1); %to improve readability below
                    
                    curSet(pntr).spikes = sp(sp>=(offset-preLag) & sp<=(handles.snipTimes(motifIdx(k),2)-preLag)) - (offset-preLag);
                    curSet(pntr).noise = noi(noi>=(offset-preLag) & noi<=(handles.snipTimes(motifIdx(k),2)-preLag)) - (offset-preLag);
                    
                    %Fuck-up check
                    if max(curSet(pntr).spikes) > (handles.snipTimes(motifIdx(k),2)-handles.snipTimes(motifIdx(k),1))
                        display('Something is fucked up... Spike times excede motif length.')
                        beep
                        return
                    end
                    
                    %Only advance pointer if the spikes file isn't empty
                    if ~isempty(curSet(pntr).spikes)
                        pntr = pntr + 1;
                    else
                        display(['No spikes found for file ' curSet(pntr).filename ' at ' num2str(curSet(pntr).snipTimes(1)) ' seconds']);
                    end
                end
            end
        else
            %If there is no matching number, throw an error/message and continue on
            display(['Yo! I''m missing a matching motif record for cell file ' num2str(cFilenums(j)) '. Check it out!'])
        end
        
    end
    
end

%Push the local variable to the output variable
syncData = curSet;

function handles = updateSelectors(handles)
%This function parses the syncData structure and generates the limits for the selection controls

%Update the record selection criteria
%.........the options for sorting by drug status
drugTemp = getFieldVectorCell(handles.syncData, 'drugStatus');
contents = ['All' sort(unique(drugTemp))];
set(handles.popup_drugstatus,'String',contents);
set(handles.popup_drugstatus,'Value',1);
handles.filter_DrugStat = contents(1);

%.........the options for sorting by directed status
directTemp = getFieldVectorCell(handles.syncData, 'directStatus');
contents = ['All' sort(unique(directTemp))];
set(handles.popup_directstatus,'String',contents);
set(handles.popup_directstatus,'Value',1);
handles.filter_DirectStat = contents(1);

%.........the upper and lower bounds of the filenum selector
filenumTemp = getFieldVector(handles.syncData, 'filenum');
set(handles.edit_fromRec,'String',num2str(min(filenumTemp)));
handles.filter_fromRec = min(filenumTemp);
set(handles.edit_toRec,'String',num2str(max(filenumTemp)));
handles.filter_toRec = max(filenumTemp);

%.........the upper and lower bounds of the time selector
datenumTemp = getFieldVector(handles.syncData, 'datenum');
dnumTemp = 24*(datenumTemp-floor(datenumTemp));
handles.filter_fromTime = min(dnumTemp);
handles.filter_toTime = max(dnumTemp);
set(handles.edit_startTime,'String',frac2str(handles.filter_fromTime));
set(handles.edit_endTime,'String',frac2str(handles.filter_toTime));

function handles = updateCellList(handles)
%Search the syncData for the cellnames that have matching motifs in the dataset
handles.cellList = sort(unique(getFieldVectorCell(handles.syncData, 'cellname')));

%Set font color using HTML formatting
pre = '<HTML><FONT color="';
post = '</FONT></HTML>';
for i = 1:numel(handles.cellList)
   str = [pre rgb2Hex(handles.colors(i,:)) '">' handles.cellList{i} post];
   listboxStr{i} = str;
end

set(handles.list_cellSelect,'String', listboxStr)

function handles = makePsth(handles)
%This function handles the plotting and formatting of PSTH for a given set of already aligned (or not) spikes. Ready-to-print
%spike times for the desired cells are passed in handles.alignedSpikes.

%Initialize variables
binCounts = [];
handles.binsFR = [];
handles.barsFR = [];

%Splines parameters
% bp = defaultParams;
% bp.prior_id = 'POISSON';
% bp.samp_iter = 1000;
% bp.use_logspline = 0;

%Set PSTH bins
binSize = str2double(get(handles.edit_psthBin, 'String')); %3ms by default
edges = 0:binSize:(size(handles.template,3));

%Determine which cells are selected for analysis
indx = get(handles.list_cellSelect,'Value');
selectCells = handles.cellList(indx);

%Select and clear raster axis
axes(handles.axes_stats); cla

%Calculate and plot PSTH per cell
cellId = getFieldVectorCell(handles.filtSync, 'cellname');
for i = 1:numel(selectCells)
    %Generate per-cell logical mask
    mask = strcmp(selectCells{i}, cellId);
    
    %Pool spike times (for each cell) across motifd
    sumSpikeTimes = cell2mat(handles.alignedSpikes(mask));
    
    %Bin spikes
    [binCounts(i,:),~] = histcounts(sumSpikeTimes,edges);
    
    %Convert units to mean firing rate
    numCellRend = numel(handles.alignedSpikes(mask));
    handles.binsFR(i,:) = binCounts(i,:)./(numCellRend*binSize/1000);
    
    %Smooth using BARS http://www.cnbc.cmu.edu/~rkelly/code.html
    f = [];
    plot((edges(1:end-1)+binSize/2), handles.binsFR(i,:), 'Color', handles.colors(indx(i),:)); hold on
end

%Format plot
axis tight
xlim([0, size(handles.template,3)])
ylabel('FR (Hz)')
set(gca,'Box', 'off', 'TickDir', 'out'); % 'XTickLabels', [], 'YTick', [])

function [isiX, isiY] = makeISI(handles)
%This function handles the calculation of interspike intervals for a given set of already aligned (or not) spikes.
%Plotting is handled by the calling function. Ready-to-print spike times for the desired cells are passed in
%handles.alignedSpikes. 

%Initialize variables
isiY = [];

%Set ISI pdf constants
minval = 0;
maxval = .25;
binSize = 0.0005; %bin size in (s)

%Determine which cells are selected for analysis
indx = get(handles.list_cellSelect,'Value');
selectCells = handles.cellList(indx);

%Calculate and plot PSTH per cell
cellId = getFieldVectorCell(handles.filtSync, 'cellname');
for i = 1:numel(selectCells)
    %Generate per-cell logical mask
    mask = strcmp(selectCells{i}, cellId);
    
    %Pool interspike intervals (i.e., the log difference between spiketimes) across motifs
    s = handles.alignedSpikes(mask);
    t = cellfun(@(x) (diff(x./1000)), s, 'UniformOutput', 0);
    sumISIs = cell2mat(t);
    
    %Calculate ISI distribution
    [isiX,isiY(i,:)] = epdf_cbins(sumISIs,binSize,minval,maxval);
end

function FR = makeFR(handles)
%This function handles the calculation of burst fraction for a given set of already aligned (or not) spikes.
%Plotting is handled by the calling function. Ready-to-print spike times for the desired cells are passed in
%handles.alignedSpikes. 

%Initialize variables
FR = [];

%Determine which cells are selected for analysis
indx = get(handles.list_cellSelect,'Value');
selectCells = handles.cellList(indx);

%Calculate and plot PSTH per cell
cellId = getFieldVectorCell(handles.filtSync, 'cellname');
for i = 1:numel(selectCells)
    %Generate per-cell logical mask
    mask = strcmp(selectCells{i}, cellId);
    
    %Calculate number of spikes across motifs
    s = handles.alignedSpikes(mask);
    t = cell2mat(cellfun(@(x) length(x), s, 'UniformOutput', 0));
    
    %Calculate mean firing rate (in Hz) by dividing by recording time
    FR{i} = t./(size(handles.template,3)/1000);
end

function BF = makeBurstFract(handles)
%This function handles the calculation of burst fraction for a given set of already aligned (or not) spikes.
%Plotting is handled by the calling function. Ready-to-print spike times for the desired cells are passed in
%handles.alignedSpikes. 

%Initialize variables
BF = [];

%Set burst constants
burstThresh = 150; %in Hz
thresh = 1/burstThresh; %convert to seconds

%Determine which cells are selected for analysis
indx = get(handles.list_cellSelect,'Value');
selectCells = handles.cellList(indx);

%Calculate per cell
cellId = getFieldVectorCell(handles.filtSync, 'cellname');
for i = 1:numel(selectCells)
    %Generate per-cell logical mask
    mask = strcmp(selectCells{i}, cellId);
    
    %Calculate interspike intervals (i.e., the log difference between spiketimes) across motifs
    s = handles.alignedSpikes(mask);
    t = cellfun(@(x) diff(x./1000), s, 'UniformOutput', 0);
    
    for j = 1:numel(t)
        %Get indices of spikes in burst
        HF_index = find(t{j}<thresh);
        spikes_adj = length(find(diff(HF_index)==1));
        
        %Calculate the total number of spikes in bursts
        HF_spikes = 2*length(HF_index)-spikes_adj;
        
        %Calculate the burst fraction for the motif
        BF{i}(j) = HF_spikes/(length(s{j}));
    end
end

function sparse = makeSparse(handles)
%This function handles the calculation of sparseness for a given set of already aligned (or not) spikes.
%Plotting is handled by the calling function. Ready-to-print spike times for the desired cells are passed in
%handles.alignedSpikes. 

%Initialize variables
binCounts = [];
binsFR = [];

%Set PSTH bins
binSize = str2double(get(handles.edit_psthBin, 'String')); %3ms by default
edges = 0:binSize:(size(handles.template,3));

%Determine which cells are selected for analysis
indx = get(handles.list_cellSelect,'Value');
selectCells = handles.cellList(indx);

%Calculate and plot PSTH per cell
cellId = getFieldVectorCell(handles.filtSync, 'cellname');
for i = 1:numel(selectCells)
    %Generate per-cell logical mask
    mask = strcmp(selectCells{i}, cellId);
    
    %Pool spike times (for each cell) across motifd
    sumSpikeTimes = cell2mat(handles.alignedSpikes(mask));
    
    %Bin spikes
    [binCounts(i,:),~] = histcounts(sumSpikeTimes,edges);
    
    %Convert units to mean firing rate
    numCellRend = numel(handles.alignedSpikes(mask));
    binsFR(i,:) = binCounts(i,:)./(numCellRend*binSize/1000);
    
    %Convert to Spiking PDF by normalizing
    pdfFR = binsFR(i,:)./sum(binsFR(i,:));
    L = log(pdfFR); L(isinf(L)) = 0; %Correct for log(0)
    sparse(i) = 1 + sum(pdfFR.*L)/log(length(pdfFR));

end

function spikeCorr = makeCorr(handles)
%This function handles the calculation of precision/correlation for a given set of already aligned (or not) spikes.
%Plotting is handled by the calling function. Ready-to-print spike times for the desired cells are passed in
%handles.alignedSpikes.

%Initialize variables
spikeCorr = [];

%Set correlation constants
gaussWidth=0.008; %in seconds
motifLength = size(handles.template,3)/1000; %in seconds
binSize = 1/handles.fs;
trainLength = floor((motifLength + eps)/ binSize + 1);

%Generate gaussian
sigma = gaussWidth / sqrt(2);
x = [gaussWidth*-4:1/handles.fs:gaussWidth*4];
gauss = (1/sqrt(2*pi)*sigma)*exp(-x.^2/(2*sigma^2));

%Determine which cells are selected for analysis
indx = get(handles.list_cellSelect,'Value');
selectCells = handles.cellList(indx);

%Calculate per cell
cellId = getFieldVectorCell(handles.filtSync, 'cellname');
for i = 1:numel(selectCells)
    %Generate per-cell logical mask
    mask = strcmp(selectCells{i}, cellId);
    
    %Retrieve the subset of spiketimes
    s = handles.alignedSpikes(mask);
    
    %Convert spike times to a matrix of spike train signals (cov w/ gaussian)
    sigs = cellfun(@(x) makeTrainSigs(x, gauss, trainLength, binSize), s, 'UniformOutput', 0);
    sigMat = cell2mat(sigs);
    
    %Calculate the Correlation matrix
    corrMat = corrcoef(sigMat);
    
    %Report only the mean correlation (minus the autocorr, i.e. the diag)
    sum_corr = sum(sum(corrMat))-sum(diag(corrMat));
    [m,n] = size(corrMat);
    spikeCorr(i) = sum_corr/((m*n)-m);

end

function sig = makeTrainSigs(x, gauss, trainLength, binSize)
%This function converts an array of spiketimes (relative to a motif start) and converts it into an analog signal. Mainly use
%for calulating pairwise correlations

%Create a binary spike train
spikeNdx = floor((x/1000 + eps) / binSize + 1);
binSpikes = zeros(trainLength, 1);
binSpikes(spikeNdx) = 1;

%Convolve binary train with the passed gaussian
sig  = conv(binSpikes, gauss);
sig = sig(floor(length(gauss)/2):end - floor((length(gauss)+1)/2)); %trim convolution artifact

function filter_Ind = filterRecords(handles)
%The purpose of this function is to take in the complete handles.syncData (which contains all the files that have matching
%records in the dataset and cell files) and return a logical index that has been appropriately parsed by the selections made
%by the user. The output is 1 if the corresponding record is to be included in the alignment; 0 if not.

%Create index for drug status
drugTemp = getFieldVectorCell(handles.syncData, 'drugStatus');
if strcmp(handles.filter_DrugStat, 'All')
    drugInd = true(length(drugTemp), 1)';
else
    drugInd = strcmp(handles.filter_DrugStat, drugTemp);
end

%Create index for directed status
directTemp = getFieldVectorCell(handles.syncData, 'directStatus');
if strcmp(handles.filter_DirectStat, 'All')
    directInd = true(length(directTemp), 1)';
else
    directInd = strcmp(handles.filter_DirectStat, directTemp);
end

%Create index for the file number limits
filenumTemp = getFieldVector(handles.syncData, 'filenum');
recInd = (filenumTemp>=handles.filter_fromRec & filenumTemp<=handles.filter_toRec);

%Create index for the filetime limits
datenumTemp = getFieldVector(handles.syncData, 'datenum');
dnumTemp = 24*(datenumTemp-floor(datenumTemp));
timeInd = (dnumTemp>=handles.filter_fromTime & dnumTemp<=handles.filter_toTime);

%Combine all indices to master filter index
filter_Ind = bitand(bitand(bitand(drugInd, directInd), recInd), timeInd);

function cellfilter_Ind = filterCells(handles)
%The purpose of this function is to take in the complete handles.syncData (which contains all the files that have matching
%records in the dataset and cell files) and return a logical index that has been appropriately parsed by the cell selections
%made by the user. The output is 1 if the corresponding record is to be included in the alignment; 0 if not.

%Retrieve the user selections from the GUI
selectCells = handles.cellList(get(handles.list_cellSelect, 'Value'));

%Extract the cellnames from the SyncData
cellId = getFieldVectorCell(handles.syncData, 'cellname');

%Index is true if the cellId has a match in the selections
cellfilter_Ind = ismember(cellId, selectCells);

function handles = showSpikes(handles)
%This function handles the plotting and formatting of already aligned (or not) spikes to the rasters axis. Ready-to-print
%spike times for the desired cells are passed in handles.alignedSpikes.

%Initialize variables
numRows = length(handles.alignedSpikes);
lineSize = 0.8;

%Generate index to plot cells by color
cellId = getFieldVectorCell(handles.filtSync, 'cellname');
colorIdx = zeros(size(cellId));
for i = 1:numel(handles.cellList)
    mask = strcmp(handles.cellList{i}, cellId);
    colorIdx(mask) = i;
end

%Sort data by user selections (only affects display)
if get(handles.check_sortdirect, 'Value')
    %Sort directed/undirected
    cellTime = getFieldVectorCell(handles.filtSync, 'directStatus');
    [~, idx] = sort(cellTime);
    
elseif get(handles.check_sortdrugs, 'Value')
    %Sort by drugs
    cellTime = getFieldVectorCell(handles.filtSync, 'drugStatus');
    [~, idx] = sort(cellTime);
    
else
    %Chronological
    cellTime = getFieldVector(handles.filtSync, 'datenum');
    [~, idx] = sort(cellTime);
end

%Select and clear raster axis
axes(handles.axes_raster); cla

%Plot each row of spikes as a single line object (much faster)
for i = 1:numRows
    [tSx, tSy]  = rasterLine(handles.alignedSpikes{idx(i)}', (i - (lineSize/2)), (i + (lineSize/2)));
    line(tSx, tSy, 'color', handles.colors(colorIdx(idx(i)),:)); hold on
end
handles.idx = idx;

%Format plot
axis tight; axis ij
xlim([0, size(handles.template,3)])
set(gca,'Box', 'off', 'TickDir', 'out', 'XAxisLocation', 'top', 'XTickLabels', [], 'YTick', [])

%Set callback for spotlight
set(gca, 'ButtonDownFcn', @cb_click);

function handles = figplotSpikes(handles)
%This function handles the plotting and formatting of already aligned (or not) spikes to the rasters axis. Ready-to-print
%spike times for the desired cells are passed in handles.alignedSpikes.

%Initialize variables
numRows = length(handles.alignedSpikes);
lineSize = 0.8;
numRends = 10;

%Generate index to plot cells by color
cellId = getFieldVectorCell(handles.filtSync, 'cellname');
colorIdx = zeros(size(cellId));
for i = 1:numel(handles.cellList)
    mask = strcmp(handles.cellList{i}, cellId);
    colorIdx(mask) = i;
end

if numRends > 0
    %Generate a list of the renditons to plot
    q = unique(colorIdx);
    for i = 1:numel(q)
        cellEdges(i) = min(find(colorIdx == q(i), 1, 'first'));
    end
    cellEdges(end+1) = numel(colorIdx) + 1;
    
    cellVec = [];
    for i = 1:numel(q)
        cellVec = [cellVec, cellEdges(i):min([(cellEdges(i)+numRends), (cellEdges(i+1)-1)])];
    end
end

%Sort data by user selections (only affects display)
if get(handles.check_sortdirect, 'Value')
    %Sort directed/undirected
    cellTime = getFieldVectorCell(handles.filtSync, 'directStatus');
    [~, idx] = sort(cellTime);
    
elseif get(handles.check_sortdrugs, 'Value')
    %Sort by drugs
    cellTime = getFieldVectorCell(handles.filtSync, 'drugStatus');
    [~, idx] = sort(cellTime);
    
else
    %Chronological
    cellTime = getFieldVector(handles.filtSync, 'datenum');
    [~, idx] = sort(cellTime);
end

%Select and clear raster axis
hFig = figure(404); clf

%Plot song template
subplot(100, 1, 1:18)
if (~isfield(handles,'template'))
    set(handles.text_message, 'String', 'Can''t find the template data. Try reloading datasets and try again.')
else
    %Show the spectrogram
    if ismatrix(handles.template)
        spec = handles.template;
    else
        spec = squeeze(handles.template(1,:,:));
    end
    imagesc(-1*spec); colormap(jet)
    axis xy; axis tight;
    set(gca,'Box', 'off', 'TickDir', 'out',  'XTickLabels', [], 'YTick', [])
    hold on
    xs = xlim;
end

ND = char(handles.datatitles(1));
title([ND(1:6) ' ' ND(8:13)])

%Plot syllable bars
subplot(100, 1, 20:22)

%Plot background
patch([xs(1), xs(2), xs(2), xs(1)], [0, 0, 1, 1], [0.5, 0.5, 0.5]); hold on

%Plot syllable bars
times = squeeze(handles.templatesyllBreaks(1,:,:));
for i = 1:size(times,1)
    patch([times(i,1), times(i,2), times(i,2), times(i,1)], [0, 0, 1, 1], 'k');
end
axis tight; xlim(xs), ylim([0,1])
set(gca,'Box', 'off', 'XTickLabels', [], 'YTick', [])

%Plot rasters
subplot(100, 1, 25:80)

%Plot each row of spikes as a single line object (much faster)
n = 1;
for i = 1:numRows
    %Should you plot it?
    plotit = false;
    if numRends < 0
        plotit = true;
    elseif ismember(i, cellVec)
        plotit = true;
    end
    
    %Conditional plot
    if plotit
        [tSx, tSy]  = rasterLine(handles.alignedSpikes{idx(i)}', (n - (lineSize/2)), (n + (lineSize/2)));
        line(tSx, tSy, 'color', handles.colors(colorIdx(idx(i)),:)); hold on
        n = n +1;
    end
end

%Format plot
axis tight; axis ij
xlim([0, size(handles.template,3)])
set(gca,'Box', 'off', 'TickDir', 'out', 'XTick', [], 'YTick', [])

%Plot PSTH
subplot(100, 1, 85:100)

%Set PSTH bins
binSize = str2double(get(handles.edit_psthBin, 'String')); %3ms by default
edges = 0:binSize:(size(handles.template,3));

%Determine which cells are selected for analysis
indx = get(handles.list_cellSelect,'Value');

for i = 1:numel(indx)
    %Plot to axis
    plot((edges(1:end-1)+binSize/2), handles.binsFR(i,:), 'Color', handles.colors(indx(i),:)); hold on
end

%Format plot
axis tight
xlim([0, size(handles.template,3)])
ylabel('FR (Hz)'); xlabel('Time (ms)')
set(gca,'Box', 'off', 'TickDir', 'out'); % 'XTickLabels', [], 'YTick', [])

set(hFig, 'Units', 'inches', 'Position', [2, 2, 5, 6])

function handles = updateSpotlight(handles)

record1 = handles.spotNum;

%Plot guide line to show what's been taken
axes(handles.axes_raster)
hold on
if isfield(handles, 'spotHand')
    delete(handles.spotHand')
end
xs = xlim;
handles.spotHand = plot(xs, [record1, record1],'k');
hold off

% Sort data given the chosen options
% [RemapIndx, ~] = sortData(handles);

% %row number
% sortChosen = handles.chosenStartSeq(RemapIndx,:);
% record2 = sortChosen(record1,1);    %File number
% startSyll = sortChosen(record1,2);   %Start syllable
% 
% %Find which file this corresponds to
% filt_keys = handles.keys(handles.filtInd==1);
% sel_filename = filt_keys(record2)
% startSyll
% 
% %Get the corresponding data structure
% filt_dataset = handles.dataset(handles.filtInd==1);
% sel_datafile = cell2mat(filt_dataset(record2));
% 
% audio = cell2mat(handles.data.PPaudio(record1));
% 
% %This section is just to format the filename for display as a figure title.
% splits = regexp(sel_filename,'_','split');
% pnt_lab =[];
% for i = 1:size(splits{1},2)-1
%     pnt_lab = [pnt_lab char(splits{1}(i)), '\_'];
% end
% pnt_lab = [pnt_lab char(splits{1}(i+1))];
% 
% if get(handles.check_singSpot,'Value')
%     figure(5656); clf
% else
%     figure; clf
% end
% displaySpecgramQuick(audio, handles.fs, [0,10000],[],0);
% title(['Segment ' pnt_lab ' / ' num2str(startSyll)])
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')

function handles = displayTemplate(handles)
%Displays spectrogram of alignment template at top. If alignment has been selected, the alignment anchors are indicated as
%well.
if (~isfield(handles,'template'))
    set(handles.text_message, 'String', 'Can''t find the template data. Try reloading datasets and try again.')
else
    %Show the spectrogram
    if ismatrix(handles.template)
        spec = handles.template;
    else
        spec = squeeze(handles.template(1,:,:));
    end
    axes(handles.axes_template); cla
    imagesc(-1*spec); colormap(jet)
    axis xy; axis tight;
    set(gca,'Box', 'off', 'TickDir', 'out', 'YTick', [])
    hold on
    
    %Plot the anchors, if selected
%     if handles.rastered
%         xs = xlim;
%         ys = ylim;
%     end
    
end

function handles = clearAllDataAxes(handles)
% Clears out all of the axes that show changeable data
axes(handles.axes_template)
cla;
axes(handles.axes_raster)
cla;
axes(handles.axes_stats)
cla;

function warpedOut = getWarpedStarts(path,anchors)
[m, n] = size(anchors);
warpedOut = zeros(m,n);

for i = 1:m
    for j = 1:n
        ind = find(path(:,1)==anchors(i,j));
        warpedOut(i,j) = round(mean(path(ind,2)));
    end
end

function hexStr = rgb2Hex(rgbColour)
%Convert a decimal RGB color code to HEX code
hexStr = reshape(dec2hex(floor(255*rgbColour), 2)', 1, 6);

function strTime = frac2str(hourFrac)
h = floor(hourFrac);
m = round(60*(hourFrac - h));
strTime = [num2str(h) ':' num2str(m,'%02.2u')];

function hourFrac = str2frac(strTime)
sp= regexp(strTime, ':', 'split');
hourFrac = str2double(sp{1}) + str2double(sp{2})/60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function push_addDataset_Callback(hObject, eventdata, handles)

%Get user input for where to find the aligned dataset
[file,path] = uigetfile('*.mat', 'Choose an aligned dataset to load:', 'MultiSelect', 'on');

if isequal(file,0) || isequal(path,0)
    set(handles.text_message,'String','Dataset location invalid or cancelled. Pick a valid file.');
else
    %Initialize the main vars
    handles.filenames = [];
    handles.filenums = [];
    handles.p = [];
    handles.q = [];
    handles.pFull = [];
    handles.snipTimes = [];
    handles.drugStatus = [];
    handles.directStatus = [];
    handles.sequence = [];
    handles.templatesyllBreaks = [];
    handles.template = [];
    handles.AlignType = [];
    handles.dataFlag = [];
    
    %Reset GUI
    handles = clearAllDataAxes(handles);
    set(handles.list_datasets, 'String', handles.datatitles);
    
    %Handle single vs multiple datasets
    if iscell(file) %proxy for multiple files
        %Sort datasets alphabeticly to maintain chronology
        file = sort(file);
        numFiles = size(file,2);
    else
        %Convert string to cell
        t{1} = file; file = t;
        numFiles = 1;
    end
    
    %Cycle through each selected file and extract the needed info
    for i = 1:numFiles
        %Keep a record of where this stuff came from
        handles.datasetFilename{i} = [path,file{i}];
        handles.datatitles{i} = file{i};
    
        %Load the dataset from file (Save time and don't import the audio)
        load(handles.datasetFilename{i}, 'sequence', 'filenames', 'filenums', 'data');
        
        %Shuffle loaded variables into the handles structure (filed per-dataset)
        sp = regexp(sequence,' ', 'split');
        handles.sequence(i,:) = str2num(sp{1});
        handles.templatesyllBreaks(i,:,:) = data.templatesyllBreaks;
        handles.template(i,:,:) = data.template;
        handles.AlignType{i} = data.AlignType;
        
        %Parse out the identifying fields of the loaded variables (filed per-dataset)
        handles.filenames = [handles.filenames; filenames'];
        handles.filenums = [handles.filenums; filenums'];
        
        %Parse out the path fields of the loaded variables (filed per-dataset)
        handles.p = [handles.p; data.p'];
        handles.q = [handles.q; data.q'];
        handles.pFull = [handles.pFull; data.pFull'];
        handles.snipTimes = [handles.snipTimes; data.snipTimes];
        handles.dataFlag = [handles.dataFlag; data.flag];
        
        %Parse out the sortable fields of the loaded variables
        handles.drugStatus = [handles.drugStatus; data.drugsStatus'];
        handles.directStatus = [handles.directStatus; data.directStatus'];
        
        %Delete the loaded variables to prevent accidental reuse
        clear('sequence', 'filenames', 'filenums', 'data')
    end

    %Verify that all datasets were aligned to the same template; if not dump to cmd
    check = sum(sum(std(handles.templatesyllBreaks,0,1)));
    if check ~= 0
        set(handles.text_message, 'String', 'Are you sure that the datasets are aligned to the same template? Something seems off... the sylBreaks don''t match.')
        return
    end
    
    %Set the message to show annotation replaced
    set(handles.list_datasets, 'String', handles.datatitles);
    set(handles.text_message, 'String', 'New datasets added to the list.');
end

    %Change indicator/button color
    if ~isempty(handles.datatitles) && ~isempty(handles.celltitles)
        set(handles.push_syncData,'BackgroundColor', [46,208,40]./255) %green-ish
    end

guidata(hObject, handles);

function push_clearDataset_Callback(hObject, eventdata, handles)
button = questdlg('Are you sure you want to clear all datsets?','Clear datasets?','Clear''em All!','Nooooo!','Nooooo!');

if (strcmp(button,'Clear''em All!'))
    %Clear out locations
    handles.datasetFilename = [];
    handles.datatitles = [];
    
    %Clear the main structure of the cell file
    handles.filenames = [];
    handles.filenums = [];
    handles.p = [];
    handles.q = [];
    handles.pFull = [];
    handles.snipTimes = [];
    handles.drugstatus = [];
    handles.directstatus = [];
    handles.sequence = [];
    handles.templatesyllBreaks = [];
    handles.template = [];
    handles.AlignType = [];

    %Set the message to show annotation replaced
    set(handles.list_datasets, 'String', handles.datatitles);
    handles = clearAllDataAxes(handles);
    
    set(handles.text_message,'String','All datasets cleared.');
    
    %Change indicator/button color
    set(handles.push_syncData,'BackgroundColor', [236, 47, 51]./255) %red-ish
else
    set(handles.text_message,'String','Datasets are unchanged.');
end

guidata(hObject, handles);

function push_addCell_Callback(hObject, eventdata, handles)
%Get user input for where to find the cell files
[file,path] = uigetfile('*.mat', 'Choose cell files to load:', 'MultiSelect', 'on');

if isequal(file,0) || isequal(path,0)
    set(handles.text_message,'String','Location(s) invalid or cancelled. Pick a valid file.');
else
    
    %Initialize the main vars
    handles.cell = [];
    
    %Reset GUI
    handles = clearAllDataAxes(handles);
    set(handles.list_datasets, 'String', handles.datatitles);
    
    %Handle single vs multiple datasets
    if iscell(file) %proxy for multiple files
        %Sort datasets alphabeticly to maintain chronology
        file = sort(file);
        numFiles = size(file,2);
    else
        %Convert string to cell
        t{1} = file; file = t;
        numFiles = 1;
    end
    
    %Cycle through each selected file and extract the needed info
    for i = 1:numFiles
        %Keep a record of where this stuff came from
        handles.cellFilename{i} = [path,file{i}];
        handles.celltitles{i} = file{i};
    
        %Load the dataset from file (Save time and don't import the audio)
        load(handles.cellFilename{i}, 'presentCell');
        
        %Shuffle loaded variables into the handles structure
        handles.cell{i} = presentCell;
        
        %Delete the loaded variables to prevent accidental reuse
        clear('presentCell')
    end

    %Set the message to show annotation replaced
    set(handles.list_cells, 'String', handles.celltitles);
    set(handles.text_message, 'String', 'New cells added to the list.');
    
    %Change indicator/button color
    if ~isempty(handles.datatitles) && ~isempty(handles.celltitles)
        set(handles.push_syncData,'BackgroundColor', [46,208,40]./255) %green-ish
    end
    
end

guidata(hObject, handles);

function push_clearCells_Callback(hObject, eventdata, handles)
button = questdlg('Are you sure you want to clear all cell files?','Clear cells?','Clear''em All!','Nooooo!','Nooooo!');

if (strcmp(button,'Clear''em All!'))
    %Clear out locations
    handles.cellFilename = [];
    handles.celltitles = [];

    %Clear the main structure of the cell file
    handles.cell = [];

    %Set the message to show annotation replaced
    set(handles.list_cells, 'String', handles.celltitles)
    handles = clearAllDataAxes(handles);
    set(handles.text_message,'String','All cell files cleared.');
    
    %Change indicator/button color
    set(handles.push_syncData,'BackgroundColor', [236, 47, 51]./255) %red-ish
else
    set(handles.text_message,'String','Cell files are unchanged.');
end

guidata(hObject, handles);

function push_syncData_Callback(hObject, eventdata, handles)
%If a datasets and cell files are already set, cross-reference the two to generate a lookup table and populate the record selection tools
if ~isempty(handles.datatitles) && ~isempty(handles.celltitles)
    %Initialize variables
    
    
    %Generate lookup table to match motifs to spikes
    handles.syncData = crossRef2(handles);
    
    %Update Selection controls
    handles = updateSelectors(handles);
    
    %Update cell list
    handles = updateCellList(handles);
    
    %Plot template to the display
    handles = displayTemplate(handles);
    
    %Update status message
    set(handles.text_message,'String','Syncing complete. Get down to it.');
    
else
    set(handles.text_message,'String','Must load dataset and cell file before syncing.');
end

guidata(hObject, handles);

function list_cells_Callback(hObject, eventdata, handles)

function list_cells_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function list_datasets_Callback(hObject, eventdata, handles)

function list_datasets_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function list_cellFilenames_Callback(hObject, eventdata, handles)

function list_cellFilenames_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function list_dataFilenames_Callback(hObject, eventdata, handles)

function list_dataFilenames_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Record Selection and Sorting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_fromRec_Callback(hObject, eventdata, handles)
%Retrieve from syncData
filenumTemp = getFieldVector(handles.syncData, 'filenum');

%Limit check the input and update
if str2double(get(handles.edit_fromRec,'String')) >= min(filenumTemp) && str2double(get(handles.edit_fromRec,'String')) <= str2double(get(handles.edit_toRec,'String'))
    handles.filter_fromRec = str2double(get(handles.edit_fromRec,'String'));
else
    set(handles.edit_fromRec,'String',num2str(min(filenumTemp)));
    handles.filter_fromRec = min(filenumTemp);
end
guidata(hObject, handles);

function edit_fromRec_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_toRec_Callback(hObject, eventdata, handles)
%Retrieve from syncData
filenumTemp = getFieldVector(handles.syncData, 'filenum');

%Limit check the input and update
if str2double(get(handles.edit_toRec,'String')) <= max(filenumTemp) && str2double(get(handles.edit_toRec,'String')) >= str2double(get(handles.edit_fromRec,'String'))
    handles.filter_toRec = str2double(get(handles.edit_toRec,'String'));
else
    set(handles.edit_toRec,'String',num2str(max(filenumTemp)));
    handles.filter_toRec = max(filenumTemp);
end
guidata(hObject, handles);

function edit_toRec_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_startTime_Callback(hObject, eventdata, handles)
%Convert str input to hour fraction
t = get(handles.edit_startTime, 'String');
startT = str2frac(t);

t = get(handles.edit_endTime, 'String');
endT= str2frac(t);

datenumTemp = getFieldVector(handles.syncData, 'datenum');
dnumTemp = 24*(datenumTemp-floor(datenumTemp));

%Limits checking
if startT >= min(dnumTemp) && startT <= endT
    handles.filter_fromTime = startT;
else
    set(handles.edit_startTime,'String',frac2str((min(dnumTemp))));
    handles.filter_fromTime = min(dnumTemp);
end

guidata(hObject, handles);

function edit_startTime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_endTime_Callback(hObject, eventdata, handles)
%Convert str input to hour fraction
t = get(handles.edit_startTime, 'String');
startT = str2frac(t);

t = get(handles.edit_endTime, 'String');
endT= str2frac(t);

datenumTemp = getFieldVector(handles.syncData, 'datenum');
dnumTemp = 24*(datenumTemp-floor(datenumTemp));

%Limits checking
if endT <= max(dnumTemp) && startT <= endT
    handles.filter_toTime = endT;
else
    set(handles.edit_endTime,'String',frac2str((max(dnumTemp))));
    handles.filter_toTime = max(dnumTemp);
end

guidata(hObject, handles);

function edit_endTime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_directstatus_Callback(hObject, eventdata, handles)

contents = get(handles.popup_directstatus,'String');
handles.filter_DirectStat = contents(get(handles.popup_directstatus,'Value'));

guidata(hObject, handles);

function popup_directstatus_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_drugstatus_Callback(hObject, eventdata, handles)

contents = get(handles.popup_drugstatus,'String');
handles.filter_DrugStat = contents(get(handles.popup_drugstatus,'Value'));

guidata(hObject, handles);

function popup_drugstatus_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function check_sortdrugs_Callback(hObject, eventdata, handles)
%Don't allow uncheck function; must select another checkbox to unselect
if ~get(handles.check_sortdrugs,'Value')
    set(handles.check_sortdrugs,'Value',1);
end

%Reset the linked checkboxes
set(handles.check_sortdirect,'Value',0);
set(handles.check_sortchrono,'Value',0);

%Replot
handles = showSpikes(handles);

guidata(hObject, handles);

function check_sortdirect_Callback(hObject, eventdata, handles)
%Don't allow uncheck function; must select another checkbox to unselect
if ~get(handles.check_sortdirect,'Value')
    set(handles.check_sortdirect,'Value',1);
end

%Reset the linked checkboxes
set(handles.check_sortdrugs,'Value',0);
set(handles.check_sortchrono,'Value',0);

%Replot
handles = showSpikes(handles);

guidata(hObject, handles);

function check_sortchrono_Callback(hObject, eventdata, handles)
%Don't allow uncheck function; must select another checkbox to unselect
if ~get(handles.check_sortchrono,'Value')
    set(handles.check_sortchrono,'Value',1);
end

%Reset the linked checkboxes
set(handles.check_sortdrugs,'Value',0);
set(handles.check_sortdirect,'Value',0);

%Replot
handles = showSpikes(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alignment Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function push_alignIt_Callback(hObject, eventdata, handles)
%Check prereqs for doing the alignement
if ~isfield(handles, 'syncData') || isempty(handles.syncData)
    set(handles.text_message, 'String', 'You gotta sync your datasets before aligning, fool! Hit that shit and try again.')
    return
end

if isempty(get(handles.list_cellSelect,'Value'))
    set(handles.text_message, 'String', 'What the hell are we doing here? You gotta select one or more cells to align.')
    return
end

%Gather processing selections
alignType = get(handles.popup_alignType,'Value');


%Apply selection filters
recfilter_Ind = filterRecords(handles); %Condition Selection
cellfilter_Ind = filterCells(handles); %Cell Selection
handles.filtSync = handles.syncData(recfilter_Ind & cellfilter_Ind);

%Align spikes
handles.alignedSpikes = [];
if alignType == 1 %Local Linear
    for i = 1:length(handles.filtSync)
        %Create sparse local linear path from interval lengths
        syllBreaks = getWarpedStarts(handles.filtSync(i).pFull, squeeze(handles.templatesyllBreaks(1,:,:)));
        a = squeeze(handles.templatesyllBreaks(1,:,:))'; path_t = [1; a(:); handles.filtSync(i).pFull(end,1)];
%         b = syllBreaks'; path_r = [1; b(:); handles.filtSync(i).pFull(end,2)];
        b = syllBreaks'; path_r = [1; b(:); ceil(diff(handles.filtSync(i).snipTimes)*1000)];
        sp_path = [path_t,path_r];
        
        handles.alignedSpikes{i} = alignSeriesPP(1000*handles.filtSync(i).spikes,sp_path);
    end
elseif alignType == 2 %Global Linear
    for i = 1:length(handles.filtSync)
        %Create sparse global linear path from interval lengths
        syllBreaks = getWarpedStarts(handles.filtSync(i).pFull, squeeze(handles.templatesyllBreaks(1,:,:)));
        a = squeeze(handles.templatesyllBreaks(1,:,:))'; path_t = [1; a(1,1); a(end,end); handles.filtSync(i).pFull(end,1)];
%         b = syllBreaks'; path_r = [1; b(1,1); b(end,end); handles.filtSync(i).pFull(end,2)];
        b = syllBreaks'; path_r = [1; b(1,1); b(end,end); ceil(diff(handles.filtSync(i).snipTimes)*1000)];
        sp_path = [path_t,path_r];
        
        handles.alignedSpikes{i} = alignSeriesPP(1000*handles.filtSync(i).spikes,sp_path);
    end
elseif alignType == 3 %Align Start
    %Get the user's selection from the GUI
    alignPnt = 2;
    
    for i = 1:length(handles.filtSync)
        %Create sparse local linear path from interval lengths
        syllBreaks = getWarpedStarts(handles.filtSync(i).pFull, squeeze(handles.templatesyllBreaks(1,:,:)));
        a = squeeze(handles.templatesyllBreaks(1,:,:))'; path_t = [1; a(:); handles.filtSync(i).pFull(end,1)];
%         b = syllBreaks'; path_r = [1; b(:); handles.filtSync(i).pFull(end,2)];
        b = syllBreaks'; path_r = [1; b(:); ceil(diff(handles.filtSync(i).snipTimes)*1000)];
        sp_path = [path_t,path_r];
        
        %Calculate the offset for edge alignment
        offset = sp_path(alignPnt,1) - sp_path(alignPnt, 2);
        
        handles.alignedSpikes{i} = (1000*handles.filtSync(i).spikes') + offset;
        
    end
elseif alignType == 4 %Full DTW
    for i = 1:length(handles.filtSync)
        handles.alignedSpikes{i} = alignSeriesPP(1000*handles.filtSync(i).spikes,handles.filtSync(i).pFull);
    end
% elseif alignType == 5 %no alignment; plot straight to axes (mostly for testing)
%     for i = 1:length(handles.filtSync)
%         handles.alignedSpikes{i} = 1000*handles.filtSync(i).spikes';
%     end
end

%Plot to GUI
handles = showSpikes(handles);

%Calculate and display PSTH
handles = makePsth(handles);

guidata(hObject, handles);

function popup_alignType_Callback(hObject, eventdata, handles)

function popup_alignType_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function check_singleSpot_Callback(hObject, eventdata, handles)

function edit_lag_Callback(hObject, eventdata, handles)

function edit_lag_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function list_cellSelect_Callback(hObject, eventdata, handles)

function list_cellSelect_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_psthBin_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function push_spikeStats_Callback(hObject, eventdata, handles)
%Do some sort of pre-req checking
if ~isfield(handles, 'alignedSpikes')
    set(handles.text_message,'String', 'You must have a set of spikes alligned before calculating statistics.')
else
    %Determine which cells are selected for analysis
    indx = get(handles.list_cellSelect,'Value');
%     selectCells = handles.cellList(indx);
    
    %Calculate Burst Fraction
    BF = makeBurstFract(handles);
   
    %Calculate Firing Rate
    FR = makeFR(handles);
    
    %Calculate Sparseness
    sparse = makeSparse(handles);
    
    %Calculate Correlation
    spikeCorr = makeCorr(handles);

    mStats = [cell2mat(cellfun(@(x) mean(x), BF, 'UniformOutput', 0))', cell2mat(cellfun(@(x) mean(x./100), FR, 'UniformOutput', 0))', sparse', spikeCorr'];
    sStats = [cell2mat(cellfun(@(x) std(x), BF, 'UniformOutput', 0))', cell2mat(cellfun(@(x) std(x./100), FR, 'UniformOutput', 0))', zeros(size(sparse')), zeros(size(spikeCorr'))];
    
    %Plot to the selected location
    if ~get(handles.check_newFig, 'Value')
        axes(handles.axes_aux); cla
    else
        figure;
    end
    
    %Plot stats as bargraphs
    b = bar(1:size(mStats,2), mStats');
    for i = 1:numel(indx)
        b(i).FaceColor = handles.colors(indx(i),:);
    end

    %Format plot
    xlim([0.5, size(mStats,2)+0.5])
    set(gca,'Box', 'off', 'TickDir', 'out', 'XTick', 1:size(mStats,2), 'XTickLabels', [{'Burst'}, {'FR/100'}, {'Sparse'}, {'Corr'}])
end

guidata(hObject, handles);

function push_corrMatrix_Callback(hObject, eventdata, handles)
%This function handles the calculation of precision/correlation matrix for a given set of already aligned (or not) spikes.
%Ready-to-print spike times for the desired cells are passed in handles.alignedSpikes.

%Do some sort of pre-req checking
if ~isfield(handles, 'alignedSpikes')
    set(handles.text_message,'String', 'You must have a set of spikes alligned before calculating statistics.')
else
    
    %Set correlation constants
    gaussWidth=0.008; %in seconds
    motifLength = size(handles.template,3)/1000; %in seconds
    binSize = 1/handles.fs;
    trainLength = floor((motifLength + eps)/ binSize + 1);
    
    %Generate gaussian
    sigma = gaussWidth / sqrt(2);
    x = [gaussWidth*-4:1/handles.fs:gaussWidth*4];
    gauss = (1/sqrt(2*pi)*sigma)*exp(-x.^2/(2*sigma^2));
    
    %Determine which cells are selected for analysis
    indx = get(handles.list_cellSelect,'Value');
    selectCells = handles.cellList(indx);
    cellId = getFieldVectorCell(handles.filtSync, 'cellname');
    
    %Build up the spiketrains to correlate
    mask = [];
    for i = 1:numel(selectCells)
        %Generate per-cell logical mask
        if isempty(mask)
            mask = strcmp(selectCells{i}, cellId);
        else
            mask = mask | strcmp(selectCells{i}, cellId);
        end
    end
    
    %Sort data by user selections (only affects display)
    if get(handles.check_sortdirect, 'Value')
        %Sort directed/undirected
        cellTime = getFieldVectorCell(handles.filtSync(mask), 'directStatus');
        [~, idx] = sort(cellTime);
        
    elseif get(handles.check_sortdrugs, 'Value')
        %Sort by drugs
        cellTime = getFieldVectorCell(handles.filtSync(mask), 'drugStatus');
        [~, idx] = sort(cellTime);
        
    else
        %Chronological
        cellTime = getFieldVector(handles.filtSync(mask), 'datenum');
        [~, idx] = sort(cellTime);
    end
    
    %Retrieve the subset of spiketimes
    s = handles.alignedSpikes(mask);
    
    %Shuffle to match displayed alignment sorting order
    s = s(idx);
    
    %Convert spike times to a matrix of spike train signals (cov w/ gaussian)
    sigs = cellfun(@(x) makeTrainSigs(x, gauss, trainLength, binSize), s, 'UniformOutput', 0);
    sigMat = cell2mat(sigs);
    
    %Calculate the Correlation matrix
    corrMat = corrcoef(sigMat);
    
    %Plot to the selected location
    if ~get(handles.check_newFig, 'Value')
        axes(handles.axes_aux); cla
    else
        figure;
    end
    
    %Plot stats as bargraphs
    imagesc(corrMat); colormap(jet)
    
    %Format plot
    axis tight
    set(gca,'Box', 'off', 'TickDir', 'out')
    
end

guidata(hObject, handles);

function push_isiDist_Callback(hObject, eventdata, handles)
%Do some sort of pre-req checking
if ~isfield(handles, 'alignedSpikes')
    set(handles.text_message,'String', 'You must have a set of spikes alligned before calculating statistics.')
else
    %Determine which cells are selected for analysis
    indx = get(handles.list_cellSelect,'Value');
    selectCells = handles.cellList(indx);
    
    %Calculate and display ISI Distributions
    [handles.isiX, handles.isiY] = makeISI(handles);
    
    %Plot ISI to the selected location
    if ~get(handles.check_newFig, 'Value')
        axes(handles.axes_aux); cla
    else
        figure;
    end
    for i = 1:numel(selectCells)
        %Plot on log scale
        semilogx(log10(handles.isiX), handles.isiY(i,:), 'Color', handles.colors(indx(i),:)); hold on
        %         plot(log10(handles.isiX), handles.isiY(i,:), 'Color', handles.colors(indx(i),:)); hold on
    end
    hold off
    
    %Format plot
    xlim([-4, 0])
    xlabel('ISI log(s)')
    ylabel('P(ISI)')
    set(gca,'Box', 'off', 'TickDir', 'out', 'XTick', [-4:0])
    
end

guidata(hObject, handles);

function push_test4_Callback(hObject, eventdata, handles)

%Function to export a spike rasters figure
handles = figplotSpikes(handles);


function push_exStats_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%
% Save to output
%%%%%%%%%%%%%%%%%%%%%
output = 'C:\Users\Tim\Desktop\TestCellFile.mat';
m = exist(output, 'file');
if m ==2
    %File already exists
    clear('cellStats');
    load(output, 'cellStats')
else
    %No file yet created
    cellStats = [];
end

%Do some sort of pre-req checking
if ~isfield(handles, 'alignedSpikes')
    set(handles.text_message,'String', 'You must have a set of spikes alligned before calculating statistics.')
else
    %Determine which cells are selected for analysis
    indx = get(handles.list_cellSelect,'Value');
    selectCells = handles.cellList(indx);
    %     selectCells = handles.cellList(indx);
    
    %Calculate Burst Fraction
    BF = makeBurstFract(handles);
    
    %Calculate Firing Rate
    FR = makeFR(handles);
    
    %Calculate Sparseness
    sparse = makeSparse(handles);
    
    %Calculate Correlation
    spikeCorr = makeCorr(handles);
    
    mStats = [cell2mat(cellfun(@(x) mean(x), BF, 'UniformOutput', 0))', cell2mat(cellfun(@(x) mean(x), FR, 'UniformOutput', 0))', sparse', spikeCorr'];
    sStats = [cell2mat(cellfun(@(x) std(x), BF, 'UniformOutput', 0))', cell2mat(cellfun(@(x) std(x), FR, 'UniformOutput', 0))', zeros(size(sparse')), zeros(size(spikeCorr'))];
    
    %Calculate and display ISI Distributions
    [handles.isiX, handles.isiY] = makeISI(handles);

    %Write the data from each cell to the output file
    cellId = getFieldVectorCell(handles.filtSync, 'cellname');
    fileId = getFieldVectorCell(handles.filtSync, 'filename');
    for i = 1:numel(selectCells)
        %Generate per-cell logical mask
        mask = strcmp(selectCells{i}, cellId);
    
        %Parse cell name
        s = regexp(selectCells{i}, '_', 'split');
        
        %Build the insertion structure
        ins = [];
        
        %Bird info
        ins.name = char(s{1});
        ins.dob = '20130801'; %YYYMMDD
        ins.tutoring = 'isolate';
        
        %Cell info
        ins.cellName = char(selectCells(i));
        ins.dor = char(s(2));
        ins.filenames = fileId(mask);
        ins.sequence = handles.sequence(1);
        
        %Cell stats
        ins.FR = mStats(i,2);
        ins.BF = mStats(i,1);
        ins.sparse =mStats(i,3);
        ins.corr = mStats(i,4);
        ins.ISI = [handles.isiX; handles.isiY(i,:)];
        ins.doc = datetime;
        
        %Conditional save
        if ~isempty(cellStats)
            ls = getFieldVectorCell(cellStats, 'cellName');
            present = ismember(ins.cellName, ls);
            
            if present
                %If it's already in the database, ask about replacing it
                button = questdlg(['Cell ' ins.cellName ' is already in the cell file. Replace current stats or skip?'],'WARNING!', 'Replace','Skip','Skip');
                
                if (strcmp(button,'Replace'))
                    pntr = find(strcmp(ins.cellName, ls), 1, 'first');
                    cellStats(pntr) = ins;
                    
                elseif (strcmp(button,'Skip'))
                    %do nothing
                    display(['Cell ' ins.cellName ' skipped.'])
                    
                end
                
            else
                %Otherwise, just append it to the end
                cellStats(end+1) = ins;
                
            end
            
        else
            %This is the first record int he cellFile
            cellStats = ins;
        end

    end
    % Sort the data by some particular field
    sortField = 'cellName';
    names = getFieldVectorCell(cellStats, sortField);
    [A, index] = sort(names);
    cellStats = cellStats(index);
    
    
    % Save the updated data to file
    save(output, 'cellStats')
    display('done')
end


guidata(hObject, handles);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spotlight Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function push_spotDelete_Callback(hObject, eventdata, handles)
if isfield(handles, 'spotHand')
    delete(handles.spotHand')
end
guidata(hObject, handles);

function push_spotUp_Callback(hObject, eventdata, handles)
%If double click then show the record number at click  
if ~isfield(handles, 'spotNum')
    handles.spotNum = 1;
else
    %Increment spotlight pointer
    handles.spotNum = handles.spotNum -1;
    
    %Limit to the size of the dataset
    handles.spotNum = max([handles.spotNum, 1]);
end

%Update the marker and the spotlight
handles = updateSpotlight(handles);

guidata(hObject, handles);

function push_spotDown_Callback(hObject, eventdata, handles)
%If double click then show the record number at click  
if ~isfield(handles, 'spotNum')
    handles.spotNum = 1;
else
    %Increment spotlight pointer
    handles.spotNum = handles.spotNum +1;
    
    %Limit to the size of the dataset
    handles.spotNum = min([handles.spotNum, numel(handles.alignedSpikes)]);
end

%Update the marker and the spotlight
handles = updateSpotlight(handles);

guidata(hObject, handles);





function edit_psthBin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function check_newFig_Callback(hObject, eventdata, handles)


a = 1;

function check_smoothPSTH_Callback(hObject, eventdata, handles)

%Update the PTSH Smoothing
a = 1;
guidata(hObject, handles);

function check_excFlags_Callback(hObject, eventdata, handles)

