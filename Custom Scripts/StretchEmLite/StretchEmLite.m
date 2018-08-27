function varargout = StretchEmLite(varargin)
% STRETCHEMLITE M-file for StretchEmLite.fig
%      Last modified by TMO 3/31/14 at 3:20am
%      STRETCHEMLITE, by itself, creates a new STRETCHEMLITE or raises the existing
%      singleton*.
%
%      H = STRETCHEMLITE returns the handle to a new STRETCHEMLITE or the handle to
%      the existing singleton*.
%
%      STRETCHEMLITE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STRETCHEMLITE.M with the given input arguments.
%
%      STRETCHEMLITE('Property','Value',...) creates a new STRETCHEMLITE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before StretchEmLite_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to StretchEmLite_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help StretchEmLite

% Last Modified by GUIDE v2.5 28-Jun-2018 12:29:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @StretchEmLite_OpeningFcn, ...
                   'gui_OutputFcn',  @StretchEmLite_OutputFcn, ...
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

function StretchEmLite_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StretchEmLite (see VARARGIN)

% Choose default command line output for StretchEmLite
handles.output = hObject;

% load('Default MU PreProcess.mat');
% handles.PreProcess = PP;

handles.annotation = mhashtable;

%linkaxes([handles.axes_audio handles.axes_neuro handles.axes_stats],'x');
set(handles.axes_audio, 'ButtonDownFcn', @cb_audio_click);
set(handles.axes_audio,'XTick',[],'YTick',[]);
set(handles.axes_stats,'YTick',[]);
set(handles.axes_template,'YTick',[]);
set(handles.axes_audioStretch,'XTick',[],'YTick',[]);
set(handles.axes_flag,'XTick',[],'YTick',[], 'Box', 'off');

%set(handles.popupSequence,'Value',1);
set(handles.popup_alignType,'Value',1);

handles.fs=44150;
handles.dirname =[];
handles.annotationFilename =[];
handles.elements = [];
handles.totFilelist = [];
handles.filelist =[];
handles.dataset =[];
handles.comTempON = false;

% Update handles structure
guidata(hObject, handles);

function varargout = StretchEmLite_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
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

guidata(hObject, handles);

function xRefFiles = XrefAnnotData(totalFiles, keys)
% This function takes in a the sorted string list "keys" which contains all
% of the .dat filenames that are required by the annotation files and
% "totalFiles" which is a structure list of all the .dat files contained
% in the loaded data folders.  "xRefFiles" will return the subset of
% totalFiles that correspond to each row of "keys".  In the event that the
% matching file cannot be found, load a dummy record in it's place

%Preallocate
fileInd = zeros(length(keys),1);
flnm = {};

%Define Dummy Record for gaps
dummy.name = 'DUMMY';
dummy.date = 'DUMMY';
dummy.bytes = 6969;
dummy.isdir = false;
dummy.datenum = 6969;
dummy.dirname = 'DUMMY';

%Extract out the filenames of all .dat files
for i = 1:length(totalFiles)
    flnm{i} = totalFiles(i).name;   
end

for i = 1:length(keys)
    index = strmatch(keys(i),flnm);
    if length(index)>1
        xRefFiles(i) = totalFiles(index(1));
    elseif isempty(index)
        xRefFiles(i) = dummy;
    else
        xRefFiles(i) = totalFiles(index);
    end
end

function filter_Ind = filterRecords(handles)
%The purpose of this function is to take in the complete handles.filelist
%(which contains all the files that are in both the data folders and in the
%annotations) and return a new list that has been appropriately parsed by
%the selections made by the user.

%Create index for drug status
if strcmp(handles.filter_DrugStat, 'All')
    drugInd = ones(length(handles.drugstatus), 1)';
else
    drugInd = strcmp(handles.filter_DrugStat, handles.drugstatus);
end

%Create index for directed status
if strcmp(handles.filter_DirectStat, 'All')
    directInd = ones(length(handles.directstatus), 1)';
else
    directInd = strcmp(handles.filter_DirectStat, handles.directstatus);
end

%Create index for the file number limits
recInd = (handles.filenums>=handles.filter_fromRec & handles.filenums<=handles.filter_toRec);

%Create index for the file number limits
timeInd = (handles.filetimes>=handles.filter_fromTime & handles.filetimes<=handles.filter_toTime);

%Combine all indices to master filter index
filter_Ind = bitand(bitand(bitand(drugInd, directInd), recInd), timeInd);

function [allSeqNumbers startSeq]=sortSequences(handles,numberSyll)
%make a cell array with all the song sequences

[allSeq startSeq counter]=getSequences(handles,numberSyll);
if (allSeq==0)
    warndlg('There are no sequences in the specified range.');
    uiwait;
    return;
end
%go through this cell array and build up an array of sequences
%numberSyll long
for i=1:numberSyll
    multiplier(i)=10^(numberSyll-i);
end
allSeqNumbers=[];
for i=1:counter-1
    allSeqNumbers(i)=dot(allSeq(i,:),multiplier);
end

function [allSeq startSeq counter]=getSequences(handles,numberSyll)
%Gets all the sequences containing numSyll.  Note that this is entirely
%dependent on the current filter settings.

%Reset placeholders
counter=1;
song_interruptions=[];

%Restrict sequence search to those files filtered by the user
filt_elements = handles.elements(handles.filtInd == 1);

%Find all pauses greater than 100ms to determine boundaries of song motifs
for i = 1:length(filt_elements)
    for j = 1:length(filt_elements{i}.segFileStartTimes)-1
        pauses(j) = filt_elements{i}.segFileStartTimes(j+1)-filt_elements{i}.segFileEndTimes(j);
        song_interruptions=find(pauses>0.2);
    end

    if length(filt_elements{i}.segFileStartTimes)>=numberSyll
        for j = 1:length(filt_elements{i}.segFileStartTimes)-numberSyll+1
            current_seq(1:numberSyll)=filt_elements{i}.segType(j:j+numberSyll-1);
            if (isempty(find(current_seq<1 | current_seq>9)) && isempty(find(song_interruptions>j-1 & song_interruptions<j+numberSyll-1)))
                allSeq(counter,1:numberSyll)=current_seq;
                startSeq(counter,:)=[i,j];%index (index in annotation file,number of starting syllable of sequences)
                counter=counter+1;
            end
        end
    end
end

if ~exist('allSeq')
    allSeq=0;
    startSeq=0;
    counter=0;
end

function handles = clearAllDataAxes(handles)
% Clears out all of the axes that show changeable data
axes(handles.axes_audio)
cla;
axes(handles.axes_stats)
cla;
axes(handles.axes_audioStretch)
cla;
axes(handles.axes_flag)
cla;

function updateFlags(handles)
%Display the flagged files array as an image in the GUI
flags = handles.data.flag;

%Plot to matrix
axes(handles.axes_flag); cla;
imagesc(flags);

axis tight;
set(gca, 'Box', 'off', 'XTick', [], 'YTick', [])

function displayTemplate(handles,sequence,syllTimes,xdata)

%syllTimes is the average syll duration
if (~isfield(handles,'templates'))
    warndlg('Load a template and try again');
    uiwait;
else
    noDisp=0;
    fs=handles.fs;
    templateDisplay(1:floor(fs*0.1))=0; %fill up buffer of 100 ms with zeros
    segTypes=getFieldVector(handles.templates.wavs,'segType');
    for i=1:length(sequence)
            presentSyll=find(segTypes==sequence(i),1,'first');
            if isempty(presentSyll)
               warndlg('no template for the chosen syllable');
               uiwait;
               templateDisplay=0;
               noDisp=1;%don't display anything
            elseif exist('syllTimes')
                noPoints= floor((syllTimes(2*i)-syllTimes(2*i-1))*fs); %number of points in proto syllable
                resampledSyll=resample(handles.templates.wavs(1,presentSyll).wav,noPoints,length(handles.templates.wavs(1,presentSyll).wav));
                templateDisplay(end+1:end+length(resampledSyll))=resampledSyll;
                if (i~=length(sequence))    
                    pause=floor((syllTimes(2*i+1)-syllTimes(2*i))*fs);
                    templateDisplay(end+1:end+pause)=0;
                end

            else    
                templateDisplay(end+1:end+length(handles.templates.wavs(1,presentSyll).wav))=handles.templates.wavs(1,presentSyll).wav;
                templateDisplay(end+1:end+floor(fs*0.015))=0;  
            end
     
    end
    if (noDisp==1)
        axes(handles.axes_template);
        cla;
    else
        
        templateDisplay(end+1:end+floor(fs*0.1))=0;
        axes(handles.axes_template);
        if (exist('xdata'))    
            displaySpecgramQuick(templateDisplay, fs, [0,10000],[],0,xdata);
        else
            displaySpecgramQuick(templateDisplay, fs, [0,10000],[],0);
        end
        set(handles.axes_template,'XTick',[]);
        set(handles.axes_template,'YTick',[]);
        xlabel('');
        ylabel('');
    end
end  

function out = normSignal(in)
%Normalize the input signal so the RMS = 1

rms = sqrt(mean(in.^2));

out = in./rms;

function filtAudio = Prep(handles,audio)
if ~isempty(audio)
    %Get filtering paramers from the GUI
    gain = str2double(get(handles.edit_sigGain, 'String'));
    HP = str2double(get(handles.edit_sigHP, 'String'));
    LP = str2double(get(handles.edit_sigLP, 'String'));
    
    %Constants for Bandpass Audio
    HP_fNorm = HP/(44150/2);
    LP_fNorm = LP/(44150/2);
    [BP_b,BP_a] = butter(4,[HP_fNorm LP_fNorm]);

    %Normalize signal RMS?
    if get(handles.check_norm, 'Value')
        audio = cellfun(@(x) normSignal(x), audio, 'UniformOutput', 0);
    end
    
    %Apply amplification and filters
    filtAudio = cellfun(@(x) filtfilt(BP_b, BP_a, gain.*(x-mean(x))), audio, 'UniformOutput', 0);
end

function SpecCube = calcPower(data)
%Generate spectrograms for each rendition
strt = 4; stp = 94; %These bins correspond to ~300-8000Hz

% Calculate STFT features for both sounds (90% window overlap)
[~,~,~,P] = cellfun(@(x) spectrogram((x/(sqrt(mean(x.^2)))),220,220-44,512,44150), data, 'UniformOutput', 0);
SpecCube = cellfun(@(x) abs(x(strt:stp,:)), P, 'UniformOutput', 0);

function [syll, motif, audio, snipTimes] = getSubSet(handles,buffer)
%Given all of the previous filtering and selecting, this function will
%parse the data from the dataset, perform the requested alignments, and
%display it on both the audio and neuro axes.

%Parse data from the dataset into two matrices (audio and neuro)
%Filter the annotation and dataset
filt_elements = handles.elements(handles.filtInd==1);
filt_dataset = handles.dataset(handles.filtInd==1);

for i=1:size(handles.chosenStartSeq,1)
    %Grab the start and end times of each syllable in the motif 
    syll(1:2:handles.numSylls*2-1,i)=filt_elements{handles.chosenStartSeq(i,1)}.segFileStartTimes(handles.chosenStartSeq(i,2):handles.chosenStartSeq(i,2)+handles.numSylls-1);
    syll(2:2:handles.numSylls*2,i)=filt_elements{handles.chosenStartSeq(i,1)}.segFileEndTimes(handles.chosenStartSeq(i,2):handles.chosenStartSeq(i,2)+handles.numSylls-1);
    
    %Grab the start and end times of each motif
    motif(1,i)= filt_elements{handles.chosenStartSeq(i,1)}.segFileStartTimes(handles.chosenStartSeq(i,2));
    motif(2,i)= filt_elements{handles.chosenStartSeq(i,1)}.segFileEndTimes(handles.chosenStartSeq(i,2)+handles.numSylls-1);

    startT = floor(((motif(1,i)-buffer)*handles.fs));
    endT = ceil(((motif(2,i)+buffer)*handles.fs));
    
   if startT < 1
        startT = 1;
    end
   if endT > length(filt_dataset{handles.chosenStartSeq(i,1)})
       endT = length(filt_dataset{handles.chosenStartSeq(i,1)});
   end
    
    %Grab the segment of data that corresponds to each selected motif
    snipTimes(i,:) = [round(startT/handles.fs,6), round(endT/handles.fs,6)];
    audio{i} = filt_dataset{handles.chosenStartSeq(i,1)}(startT:endT);

end

function template_Intermed = createSpecAve(data,buffer,starts,edgeThresh,band,offsetLim)
%Select the average-length rendition to use as the template
renditions = size(data,2);
syllBreaks = []; motifBreak = []; ints = [];

%Calculate the total power in the rendition spectrogram
s = cellfun(@(x) sum(sum(x)), data, 'UniformOutput', 0);
specPow = mean(cell2mat(s));

for i=1:renditions
    try
        %Find the motif edges
        [syllBreaks(i,:,:), motifBreak(i, :)] = getEdges(data{i},starts(:,i),buffer,0,'robust', offsetLim);
        
    catch
        %A rough fix incase one of the files fails...
        syllBreaks(i,:,:) = round(mean(syllBreaks, 1));
        motifBreak(i, :) = round(mean(motifBreak, 1));
    end
end

%calculate motif length
motifLength = motifBreak(:, 2) - motifBreak(:, 1);

%Takes the most-average length rendition and uses it alone as initial template
[~, target] = min(abs(motifLength - mean(motifLength)));
template_Initial = data{target(1)};

%Warp all spectrograms to the initial template using Cengiz's method
[~, p] = cellfun(@(x) DTWFinal(template_Initial,x,band), data, 'UniformOutput', 0);
for i=1:renditions
     audioCube(i,:,:) = alignSeries(data{i},p{i});
end

%Take the mean across all time-freq bins
Spec = squeeze(mean(audioCube,1));
[tBreaks, ~] = getEdges(Spec,starts(:,target(1)),buffer,0,'robust', offsetLim);

%Scale Power to match sample average
% template_Intermed = Spec.*(mean(specPow)/sum(sum(Spec)));  <======= THIS IS THE ORIGINAL 

%Unwarp the derived template to the so that the length of each interval is the mean across the rendition stack
%template rendition anchors
a = tBreaks'; 
tempSyllBreaks = round([1; a(:); size(Spec,2)]);

%Mean dur across intervals
ind = 1:(size(syllBreaks,2)-1);
meanOnSet = mean(squeeze(syllBreaks(:,1,1)) - ones(renditions, 1));
meanOffSet = mean(cell2mat(cellfun(@(x) size(x,2), data, 'UniformOutput', 0)) - squeeze(syllBreaks(:,end,end))');
meanSylDur = mean(diff(syllBreaks,1,3));
meanGapDur = mean(syllBreaks(:,ind+1,1) - syllBreaks(:,ind,2));

%Reconstruct mean paths
ints(1:2:(numel(meanSylDur)*2)) = meanSylDur;
ints(2:2:((numel(meanSylDur)*2)-1)) = meanGapDur;
ints = [1, meanOnSet, ints, meanOffSet];
anchors = round(cumsum(ints)');

%Construct linear paths from the derived anchors
% sp_path = [tempSyllBreaks, anchors];
sp_path = [anchors, tempSyllBreaks];

%Warp derived template
platTemp = alignSeriesSTW(Spec,sp_path);

%Scale Power to match sample average
template_Intermed = platTemp.*(specPow/sum(sum(platTemp)));

function [syllBreaks,motifBreak] = getEdges(data,starts,buffer,edgeThresh,type, offsetLim)

Crossings.Up = [];
Crossings.Down = [];
    
%Sum across frequency bins for each rendition
specPowTS = -1*sum(data);
specPowTS = mat2gray(specPowTS);

%Rerun EM until a solution is found
iterations = 1;

while (isempty(Crossings.Up) || isempty(Crossings.Down))
    gMix = gmdistribution.fit(specPowTS',2);
    [powerM,pnt] = sort(gMix.mu);
    temp = sqrt(squeeze(gMix.Sigma));
    powerStd = temp(pnt);
    thresh = powerM(1)+edgeThresh*powerStd(1); %1SDs above the silent Gaussian mean

    %Find the crossing points to estimate the onset and offsets of
    %motifs
    ind = 2:length(specPowTS);
    Crossings.Up = find(specPowTS(ind)>thresh & specPowTS(ind-1)<thresh);
    Crossings.Down = find(specPowTS(ind)<thresh & specPowTS(ind-1)>thresh);
    
    if iterations > 5
        print(['gmdistribution.fit has been run ' num2str(iterations) ' without converging on a solution.'])
        return
    end
    iterations = iterations + 1;
end

%Use the starting positions to capture the start/stop points for all
%syllables of the rendition
syllBreaks = [];
motifBreak = [];

onStarts = 1000*(starts(1:2:end)-starts(1)+buffer);
offStarts = 1000*(starts(2:2:end)-starts(1)+buffer);

 for k = 1:length(onStarts)
     [offset, pntr] = min(abs(Crossings.Up-onStarts(k)));
      if offset > offsetLim && strcmp(type,'robust')
         syllBreaks(k,1) = onStarts(k);
      else
        syllBreaks(k,1) = Crossings.Up(pntr);
      end
     
     [offset, pntr] = min(abs(Crossings.Down-offStarts(k)));
      if offset > offsetLim && strcmp(type,'robust')
         syllBreaks(k,2) = offStarts(k);
      else
        syllBreaks(k,2) = Crossings.Down(pntr);
      end
 end
motifBreak(1) = syllBreaks(1,1); motifBreak(2) = syllBreaks(end,2);

function warpedOut = getWarpedStarts(path,anchors)
[m, n] = size(anchors);
warpedOut = zeros(m,n);

for i = 1:m
    for j = 1:n
        ind = find(path(:,1)==anchors(i,j));
        warpedOut(i,j) = round(mean(path(ind,2)));
    end
end

function [RemapInx,breakpnts,RemapDrugs,RemapDirect] = sortData(handles)
%Based on the checkbox selections, create an index that orders the motifs;
%note that within categories, ordering is chronological
rendNum = size(handles.data.audioPower,1);

%Get all the 
datasetPntr = handles.chosenStartSeq(:,1);    %file number

%Get the corresponding data structure
RemapInx = [];
breakpnts = [];

filt_drugs = handles.drugstatus(handles.filtInd==1); %Strip out drugs stats for filtered files
RemapDrugs = filt_drugs(datasetPntr); %Create vector of drug stats for rendition
filt_direct = handles.directstatus(handles.filtInd==1); %Strip out drugs stats for filtered files
RemapDirect = filt_direct(datasetPntr); %Create vector of drug stats for rendition

if get(handles.check_sortdrugs,'Value')
    %Create a sorted stack for each drug stat type
    contents = get(handles.popup_drugstatus,'String');
    for i = 2:length(contents)
        bitBang = strcmp(contents(i,:),RemapDrugs);
        stack{i-1} = find(bitBang==1);
    end
    %Collapse into a single array and note breaks in data types
    for i = 1:length(stack)
        RemapInx = [RemapInx, stack{i}];
        breakpnts(i) = stack{i}(end);
    end
elseif get(handles.check_sortdirect,'Value')
    dbitBang = strcmp('Directed',RemapDirect);
    directed = find(dbitBang==1);
    ubitBang = strcmp('Undirected',RemapDirect);
    undirected = find(ubitBang==1);
    RemapInx = [undirected, directed];
    breakpnts = length(undirected);
else
    RemapInx = 1:rendNum;
end

function handles = showAlignedData(handles,sc_audio,RemapIndx,breaks)

% renditions = size(handles.data.audio,2);

%Show data on axes
axes(handles.axes_audio);
h1=imagesc(sc_audio);
set(h1, 'ButtonDownFcn', @cb_click);
hold on
for i = 1:length(breaks)
    line([0 size(handles.data.audioPower,2)],[breaks(i) breaks(i)],'Color','w','LineWidth',1)
end
hold off

%plot basic statistics across axes_stats.  Here, I'm showing the mean
%across all trials +/- standard deviation.
axes(handles.axes_stats)
cla;
handles.data.mAudio = mean(sc_audio);
plot(mean(sc_audio));
hold on
plot(mean(sc_audio)+std(sc_audio),':');
plot(mean(sc_audio)-std(sc_audio),':');
axis tight

%NO ROTATION>>> Subtract unity vector from the linear path...no need to
%interpolate or anything fancy
rotatedY = (cell2mat(handles.data.q)-cell2mat(handles.data.p))';
handles.data.pathR = rotatedY;
handles.data.pathT = handles.data.p{end};

axes(handles.axes_audioStretch)
imagesc(rotatedY(RemapIndx,:));

%Create and plot the popup figure showing different path representations
figure
subplot(2,2,1:2)
imagesc(rotatedY(RemapIndx,:),[-30,30]);

subplot(2,2,3)
plot(cell2mat(handles.data.q))
axis tight

subplot(2,2,4)
hold on
plot(rotatedY')
axis tight

%Add the mean path and unity-path on top
handles.data.pathMean = mean(rotatedY);
plot(handles.data.p{end},handles.data.pathMean,'k', 'LineWidth',3)
plot([0,max(handles.data.p{end})],[0,0],'y','LineWidth',2)
axis tight

set(handles.axes_audio,'XTick',[],'YTick',[]);
set(handles.axes_audioStretch,'XTick',[],'YTick',[]);

function strTime = frac2str(hourFrac)
h = floor(hourFrac);
m = round(60*(hourFrac - h));
strTime = [num2str(h) ':' num2str(m)];

function hourFrac = str2frac(strTime)
sp= regexp(strTime, ':', 'split');
hourFrac = str2double(sp{1}) + str2double(sp{2})/60;

function handles = updateSpotlight(handles)

record1 = handles.spotNum;

%Plot guide line to show what's been taken
axes(handles.axes_audio)
hold on
if isfield(handles, 'spotHand')
    delete(handles.spotHand')
end
xs = xlim;
handles.spotHand = plot(xs, [record1, record1],'k');
hold off

% Sort data given the chosen options
[RemapIndx, ~] = sortData(handles);

%row number
sortChosen = handles.chosenStartSeq(RemapIndx,:);
record2 = sortChosen(record1,1);    %File number
startSyll = sortChosen(record1,2);   %Start syllable

%Find which file this corresponds to
filt_keys = handles.keys(handles.filtInd==1);
sel_filename = filt_keys(record2);
% startSyll

%Get the corresponding data structure
filt_dataset = handles.dataset(handles.filtInd==1);
sel_datafile = cell2mat(filt_dataset(record2));

% audio = cell2mat(handles.data.PPaudio(record1));
audio = -handles.data.audioSpecs_LN{record1};
tempBreaks = handles.data.templatesyllBreaks(:);
syllBreaks = getWarpedStarts(handles.data.pFull{record1}, tempBreaks);

%This section is just to format the filename for display as a figure title.
splits = regexp(sel_filename,'_','split');
pnt_lab =[];
for i = 1:size(splits{1},2)-1
    pnt_lab = [pnt_lab char(splits{1}(i)), '\_'];
end
pnt_lab = [pnt_lab char(splits{1}(i+1))];

if get(handles.check_singSpot,'Value')
    figure(5656); clf
else
    figure; clf
end
% displaySpecgramQuick(audio, handles.fs, [0,10000],[],0);
% title(['Segment ' pnt_lab ' / ' num2str(startSyll)])
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')

%More useful splotlight
subplot(2,1,1); cla
imagesc(audio); hold on
axis xy; axis tight; colormap(jet)
y = ylim;
for i = 1:numel(syllBreaks)
    line([syllBreaks(i), syllBreaks(i)], y, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.5) 
end
set(gca,'Box', 'off', 'TickDir', 'out', 'YTick',[]);
title(['Segment ' pnt_lab ' / ' num2str(startSyll)])


subplot(2,1,2); cla
imagesc(-1.*handles.data.template); hold on
axis xy; axis tight; colormap(jet)
y = ylim;
for i = 1:numel(tempBreaks)
    line([tempBreaks(i), tempBreaks(i)], y, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.5) 
end
set(gca,'Box', 'off', 'TickDir', 'out', 'YTick',[]);
title('Alignment Template')
xlabel('Time (ms)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Motif Selections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_seqStats_Callback(hObject, eventdata, handles)
%go through the annotations and whenever there is a stretch of contiguous
%sequence assign the transition to the matrix. display the figure.
numberSyll=2;
[allSeq startSeq counter]=getSequences(handles,numberSyll);
numberSyllTypes=max(max(allSeq));
seqMatrix=zeros(numberSyllTypes,numberSyllTypes);   
for i=1:counter-1
    from=allSeq(i,1);
    to=allSeq(i,2);
    seqMatrix(from,to)=seqMatrix(from,to)+1;
end
figure;
h=imagesc([1 numberSyllTypes],[1 numberSyllTypes],seqMatrix);

function popup_seqSylls_Callback(hObject, eventdata, handles)

%Returns the numerical sequence that the user has selected
chosenSeq = handles.seq(get(handles.popup_seqSylls, 'Value')); %chosen sequence as number
handles.Seq = chosenSeq;

%Finds the index of the chosen sequence as they appear in the
%handles.all_seq catalog
startIndex = find(handles.allSeq==chosenSeq);

%Gives the filtered record/annotation number and the syllable number within
%the annotation
handles.chosenStartSeq = [];
for i=1:length(startIndex)
    handles.chosenStartSeq(i,:) = handles.startSeq(startIndex(i),:); 
end

%Filter the elements array 
filt_elements = handles.elements(handles.filtInd == 1);

%Determine the selected sequence and prep it for displaying the appropriate
%template on axes_template
segType = filt_elements{handles.chosenStartSeq(1,1)}.segType;
handles.sequence=[];
handles.sequence(1:handles.numSylls) = segType(handles.chosenStartSeq(1,2):handles.chosenStartSeq(1,2)+handles.numSylls-1);

%How the hell is this going to work with the warping??? 
displayTemplate(handles,handles.sequence);

%Clear out old data and templates related to the previous selections
handles.data = [];

guidata(hObject, handles);

function popup_seqSylls_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_numSylls_Callback(hObject, eventdata, handles)

if (~isfield(handles,'filter_fromRec') || ~isfield(handles,'filter_toRec') || ~isfield(handles,'filter_DrugStat') || ~isfield(handles,'filter_DirectStat'))
    warndlg('Select the records that you want to display before continuing.');
    uiwait;
elseif isempty(handles.dataset)
    set(handles.text_message,'String','There is no dataset in memory.  Create dataset before continuing.')
    set(handles.popup_numSylls,'Value',1);
else
    %Copy out the user's selection from the popup menu
    handles.numSylls = get(handles.popup_numSylls, 'Value');

    %Clear out all of the plotted data
    handles = clearAllDataAxes(handles);
    
    %Create a record index based on user filter selections
    handles.filtInd = filterRecords(handles);
    
    %Run through the now refined set of records to extract the sequences
    [handles.allSeq handles.startSeq]=sortSequences(handles, handles.numSylls); %sort the sequences of handles.numSyll in terms of how often they appear

    %Sort and display sequences by prevalance
    histx=(1:10^handles.numSylls);
    hist_seq=hist(handles.allSeq,histx);
    
    %Choose the 50 most abundant sequences to display in the sequence pop-up window
    [sortSeq, sortIndex]=sort(hist_seq,'descend');
    numSeq=length(find(sortSeq>0));
    if isempty(numSeq)
        warndlg('No such sequences found.');
        uiwait;
        return;
    elseif (numSeq<50)
        popupNum=numSeq;
    else
        popupNum=50;
    end
    popupSeqString(1,:)=[num2str(sortIndex(1)),' (n=', ' ',num2str(sortSeq(1)),')'];
    
    handles.seq=sortIndex(1);
    if (popupNum>1)
        for i=2:popupNum %display the frequency of the various sequences)
            padLength=length(num2str(sortSeq(1)))-length(num2str(sortSeq(i)));
            if (padLength==1)
                freq(i,:)=['  ',num2str(sortSeq(i))];
            elseif (padLength==2)
                freq(i,:)=['   ',num2str(sortSeq(i))];
            elseif (padLength==3)
                freq(i,:)=['    ',num2str(sortSeq(i))];
            elseif (padLength==4)
                freq(i,:)=['     ',num2str(sortSeq(i))];
            elseif (padLength==5)
                freq(i,:)=['      ',num2str(sortSeq(i))];
            elseif (padLength==6)
                freq(i,:)=['       ',num2str(sortSeq(i))];
            else
                freq(i,:)=[' ',num2str(sortSeq(i))];
            end
            handles.seq(i)=sortIndex(i);
            popupSeqString(i,:)=[num2str(sortIndex(i)),' (n=', freq(i,:),')'];
        end
    end
        
    set(handles.popup_seqSylls,'String',popupSeqString);
    set(handles.popup_seqSylls,'Value',1);
end

guidata(hObject, handles);

function popup_numSylls_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alignment Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_alignIt_Callback(hObject, eventdata, handles)
%Gather variable from GUI for later use
buffer = str2double(get(handles.edit_buffer,'String'))/1000;
band = str2double(get(handles.edit_dtwBand,'String'));
offsetLim = str2double(get(handles.edit_offsetLims,'String'));
edgeThresh = str2double(get(handles.edit_edgeThresh,'String'));
alignType = get(handles.popup_alignType,'Value');

%Capture the subset of the data set that corresponds to the selected motifs
[handles.data.syll, handles.data.motif, handles.data.audio, handles.data.snipTimes] = getSubSet(handles,buffer);

%Filter audio
[handles.data.PPaudio] = Prep(handles,handles.data.audio);

%Create Spectrogram matrix for each motif rendition
handles.data.audioSpecs = calcPower(handles.data.PPaudio);

%Convert to log units and normalize the SpecCubes
handles.data.audioSpecs_LN = cellfun(@(x) -1*log10(x), handles.data.audioSpecs, 'UniformOutput', 0);

%Generate template for the alignment
if handles.comTempON
    %Copy the common template to the active structure
    handles.data.template = handles.comTemplate;
    handles.data.templatesyllBreaks = handles.comTemplatesyllBreaks;
    handles.data.templatemotifBreaks = handles.comTemplatemotifBreaks;
else
    %Create template for alignment
    handles.data.template = createSpecAve(handles.data.audioSpecs_LN,buffer,handles.data.syll,edgeThresh,band,offsetLim);
     
    %Send template to the template bender for checking
    templateIn.template = handles.data.template;
    templateIn.threshold = edgeThresh;
    output = TemplateBender(templateIn);
    
    %Copy out of the template bender
    handles.data.template = output.template;
    handles.data.templatesyllBreaks = output.templatesyllBreaks;
    handles.data.templatemotifBreaks = [output.templatesyllBreaks(1,1), output.templatesyllBreaks(end, end)];
    set(handles.edit_edgeThresh,'String', num2str(output.threshold));
end

%Plot the spectrogram of the template on the top bar
axes(handles.axes_template)
imagesc(-1*handles.data.template)
axis xy; axis tight;
set(gca,'Box', 'off', 'TickDir', 'out', 'YTick', [])

%Reset variables
handles.data.aligned_audioCube = [];
handles.data.syllBreaks = [];

%Align the audio data to the selected template
h = waitbar(0, 'Chill out... maybe go have a smoke. I''m aligning this shit.');
[~, timeBins] = size(handles.data.template);
for i=1:size(handles.data.audio,2);
    %Generate warp path for the audio spectrogram to the template
    [~,rendLength] = size(handles.data.audioSpecs_LN{i});
    
    if get(handles.check_floatDTW,'Value')
        [~, p] = DTWsubWeightedNormBandFast(handles.data.template,handles.data.audioSpecs_LN{i}, band);
    else
        [~, p] = DTWFinal(handles.data.template,handles.data.audioSpecs_LN{i}, band);
    end
    
    %Retain full warping path
    handles.data.pFull{i} = p;
    
    %Get the syllable start/stops based on warping path
    handles.data.syllBreaks(i,:,:) = getWarpedStarts(p,handles.data.templatesyllBreaks);
    
    %Generate the alignment path based on user selections
    if alignType == 1 %Local Linear
        %Create sparse local linear path from interval lengths
        a = handles.data.templatesyllBreaks'; path_t = [1; a(:); timeBins];
        b = squeeze(handles.data.syllBreaks(i,:,:))'; path_r = [1; b(:); rendLength];
        sp_path = [path_t,path_r];
        
        %Create full path by linear interpolation
        handles.data.p{i} = (1:timeBins)';
        handles.data.q{i} = interp1(path_t,path_r,1:timeBins,'linear')';
        
    elseif alignType == 2 %Global Linear
        %Create sparse global linear path from interval lengths
        a = handles.data.templatesyllBreaks'; path_t = [1; a(1,1); a(end,end); timeBins];
        b = squeeze(handles.data.syllBreaks(i,:,:))'; path_r = [1; b(1,1); b(end,end); rendLength];
        sp_path = [path_t,path_r];
        
        %Create full path by linear interpolation
        handles.data.p{i} = (1:timeBins)';
        handles.data.q{i} = interp1(path_t,path_r,1:timeBins,'linear')';
    elseif alignType == 3 %Align Start
        
        display('Tim hasn"t coded this function yet... Sorry!')
        return
        
    elseif alignType == 4 %Full DTW
        %Create copy over the full DTW path (for programming consistency)
        sp_path = p;
        [path, ~] = alignSeriesPath(handles.data.audioSpecs_LN{i},p);
        
        %Create full path functions
%         handles.data.p{i} = p(:,1);
%         handles.data.q{i} = p(:,2);
        handles.data.p{i} = (1:timeBins)';
        handles.data.q{i} = path';
        
        if numel(handles.data.p{i}) ~= numel(handles.data.q{i})
            t854 = 1;
        end
    end
    
    %Warp the rendition spectrogram using the sparse path
    handles.data.aligned_audioCube(i,:,:) = alignSeriesSTW(handles.data.audioSpecs_LN{i},sp_path);
    
    waitbar(i/size(handles.data.audio,2))
end
close(h)

%Create power contour plots for audio
handles.data.audioPower = -1*abs(squeeze(sum(handles.data.aligned_audioCube,2)));

% Sort data given the chosen options
[RemapIndx,breaks,handles.data.drugsStatus,handles.data.directStatus] = sortData(handles);

if get(handles.check_sortdrugs,'Value')
    handles.data.drugsRemapInx = RemapIndx;
elseif get(handles.check_sortdirect,'Value')
    handles.data.directRemapInx = RemapIndx;
elseif get(handles.check_sortchrono,'Value')
    handles.data.chronoRemapInx = RemapIndx;
end

%Generate logical flag index, initially set to false
handles.data.flag = false(size(handles.data.audio,2), 1);

%Update flag display
updateFlags(handles);

sc_audio = handles.data.audioPower(RemapIndx,:);

handles = showAlignedData(handles,sc_audio,RemapIndx,breaks);

%Popup/update path deviation figure
figure(666)
plot(diff(handles.data.pathMean))
hold on
axis tight
xlabel('Template Time (ms)')
ylabel('<--Shrink        Stretch-->')
title('Mean Path Deviation from Unity Stretch')

guidata(hObject, handles);

function popup_alignType_Callback(hObject, eventdata, handles)

function popup_alignType_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function push_comTemplate_Callback(hObject, eventdata, handles)

if handles.comTempON        %Exit Common Template Mode...
    %Clear the handles structure
    handles.comTemplate = [];
    handles.comTemplatesyllBreaks = [];
    handles.comTemplatemotifBreaks = [];
    
    %Set flags for state change
    handles.comTempON = false;
    set(handles.push_comTemplate,'BackgroundColor',[0.9412, 0.9412, 0.9412]);
    set(handles.text_message,'String','Exited Common Template Mode.');
else                        %If not, go into Common Template Mode
    %Retrieve the template file
    [fname pathname] = uigetfile('*.mat','Select *.mat file containing template information');
    load([pathname fname]);
    
    %Copy the loaded variables to the handles structure
    handles.comTemplate = template;
    handles.comTemplatesyllBreaks = templatesyllBreaks;
    handles.comTemplatemotifBreaks = templatemotifBreaks;
    
    %Set flags for state change
    handles.comTempON = true;
    set(handles.push_comTemplate,'BackgroundColor','r')
    set(handles.text_message,'String',['Entered Common Template Mode. Template loaded from ' pathname fname '.']);
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
    handles.spotNum = min([handles.spotNum, size(handles.data.audio,2)]);
end

%Update the marker and the spotlight
handles = updateSpotlight(handles);

guidata(hObject, handles);

function check_singSpot_Callback(hObject, eventdata, handles)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function listbox_annotation_Callback(hObject, eventdata, handles)

function listbox_annotation_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_data_Callback(hObject, eventdata, handles)

function listbox_data_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function push_addFolder_Callback(hObject, eventdata, handles)

%Get user input for where to find the .dat files and check location
directory_name = uigetdir('Choose location of recordings to analyze:');

if isequal(directory_name,0) || ~isstr(directory_name)
    set(handles.text_message,'String','Data location invalid or cancelled. Pick a valid folder.');
    guidata(hObject, handles);
    return
end
cd(directory_name)

%Get WAV and DAT filenames, according to user selection
filelistWAV = [];
filelistDAT = [];
filelistSTM = [];
if get(handles.check_loadWAV, 'Value')
    filelistWAV = dir('*.wav');
end
if get(handles.check_loadDAT, 'Value')
    filelistDAT = dir('*.dat');
end
if get(handles.check_loadSTM, 'Value')
    filelistSTM = dir('*.stm');
end
filelist = ([filelistWAV;filelistDAT;filelistSTM]);
for i = 1:length(filelist)
    filelistnames{i} = filelist(i).name;
end
[~, indx] = sort(filelistnames);
filelist = filelist(indx);

if isempty(filelist)
    set(handles.text_message,'String',[directory_name ' does not contain any requested files. Pick a valid folder.']);
    guidata(hObject, handles);
    return
else   
    %If they don't hit cancel, update data location...

    %This is the code for handling multiple cell files. 
    if ~isempty(handles.dirname)
        button = questdlg('Do you want to add this data location to the previously loaded one(s) or replace them?','Add a data location','Add to current','Replace current','Add to current');
        if (strcmp(button,'Replace current'))
            %Overwrite the current directory handle
            handles.dirname = directory_name;
            
            %Add directory location to the file structure
            for i = 1:length(filelist)
                filelist(i).dirname = directory_name;
            end
            handles.totFilelist = filelist;
            
            %Set the message to show dataset expanded
            set(handles.text_message,'String',['Files from from ' directory_name ' are the new dataset.']);    
        elseif (strcmp(button,'Add to current'))
            if strcmp(directory_name,handles.dirname)
                set(handles.text_message,'String',[directory_name ' is already included in the dataset.']);
                guidata(hObject, handles);
                return
            else
                %handles.dirname = {handles.dirname; directory_name};
                
                if size(handles.dirname,1)==1
                    handles.dirname = {handles.dirname;directory_name};
                else
                    handles.dirname{length(handles.dirname)+1} = directory_name;
                end
            end
            
            %Add directory location to the file structure
            for i = 1:length(filelist)
                filelist(i).dirname = directory_name;
            end
            handles.totFilelist =  [handles.totFilelist; filelist];
            
            %Set the message to show dataset expanded
            set(handles.text_message,'String',['Files from from ' directory_name ' added to the dataset.']);
        end
    %This is the code for loading the first/only folder  
    else
        %Overwrite the current directory handle
        handles.dirname = directory_name;
        
        %Add directory location to the file structure
        for i = 1:length(filelist)
            filelist(i).dirname = directory_name;
        end
        handles.totFilelist = filelist;

        %Set the message to show dataset expanded
        set(handles.text_message,'String',['Files from from ' directory_name ' are the new dataset.']);
    end
    set(handles.listbox_data,'String',handles.dirname);
end

%If the annotation file is already loaded, cross reference the data
if ~isempty(handles.elements)
    %Pare down the totFilelist set to only those files that have a record
    %in the loaded annotations; dump the results in handles.filelist
    [handles.filelist] = XrefAnnotData(handles.totFilelist, handles.keys);
    set(handles.listbox_keys,'String', handles.keys)
    for i=1:length(handles.filelist)
        names{i} = handles.filelist(i).name;
        t = getFileTime(handles.filelist(i).name);
        handles.filetimes(i) = 24*(t-floor(t)); %This provides the time as the fractional hour of the day
    end
    set(handles.listbox_filelist,'String',names)
        
    %Update the record selection criteria
    %.........the options for sorting by drug status
    contents = ['All' sort(unique(handles.drugstatus))];
    set(handles.popup_drugstatus,'String',contents);
    set(handles.popup_drugstatus,'Value',1);
    handles.filter_DrugStat = contents(1);

    %.........the options for sorting by directed status
    contents = ['All' sort(unique(handles.directstatus))];
    set(handles.popup_directstatus,'String',contents);
    set(handles.popup_directstatus,'Value',1);
    handles.filter_DirectStat = contents(1);

    %.........the upper and lower bounds of the filenum selector
    set(handles.edit_fromRec,'String',num2str(min(handles.filenums)));
    handles.filter_fromRec = min(handles.filenums);
    set(handles.edit_toRec,'String',num2str(max(handles.filenums)));
    handles.filter_toRec = max(handles.filenums);
    
    %.........the upper and lower bounds of the time selector
    handles.filter_fromTime = min(handles.filetimes);
    handles.filter_toTime = max(handles.filetimes);
    set(handles.edit_startTime,'String',frac2str(handles.filter_fromTime));
    set(handles.edit_endTime,'String',frac2str(handles.filter_toTime));
    
end

guidata(hObject, handles);

function push_addAnnot_Callback(hObject, eventdata, handles)
% hObject    handle to push_addAnnot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get user input for where to find the annotation file
[file,path] = uigetfile('*.mat', 'Choose an audio annotation .mat file to load:');

if isequal(file,0) || isequal(path,0)
    set(handles.text_message,'String','Annotation location invalid or cancelled. Pick a valid file.');
else
    %if they don't hit cancel, update annotation location and load the file
    if ~isempty(handles.annotationFilename)
        button = questdlg('Do you want to add this annotation to the current list or replace them?','Add an annotation','Add to current','Replace current','Add to current');
        if (strcmp(button,'Replace current'))
            handles.annotationFilename = [path,file];
            handles.annotitles = file;
            annotation = aaLoadHashtable(handles.annotationFilename);

            %Copy out the two main structures of the annotation file
            handles.elements = annotation.elements;
            handles.keys = annotation.keys;
            
            %Sort the keys alphabetically (which is also temporally) and
            %use the index to sort the elements vector as well
            [handles.keys, ind] = sort(handles.keys);
            handles.elements = handles.elements(ind);
            
            %Parse out the sortable fields of the elements structure
            handles.filenums=getAnnotationVector(handles.elements,'filenum');
            handles.drugstatus=getAnnotationVector(handles.elements,'drugstatus');
            handles.directstatus=getAnnotationVector(handles.elements,'directstatus');
            
            %Set the message to show annotation replaced
            set(handles.listbox_annotation,'String',handles.annotitles)
            set(handles.text_message,'String',[file ' is the new annotation.']);
            
        elseif (strcmp(button,'Add to current'))
            if strcmp(file,handles.annotitles)
                set(handles.text_message,'String',[file ' is already included in the dataset.']);
                guidata(hObject, handles);
                return
            end
            
            if size(handles.annotitles,1)==1
                handles.annotationFilename = {handles.annotationFilename; [path,file]};
                handles.annotitles = {handles.annotitles;file};
            else
                handles.annotationFilename{length(handles.annotationFilename)+1} = [path,file];
                handles.annotitles{length(handles.annotitles)+1} = file;
            end
            
            annotation = aaLoadHashtable([path,file]);
            
            %Copy out the two main structures of the annotation file
            handles.elements = [handles.elements annotation.elements];
            handles.keys = [handles.keys annotation.keys];
            
            %Sort the keys alphabetically (which is also temporally) and
            %use the index to sort the elements vector as well
            [handles.keys, ind] = sort(handles.keys);
            handles.elements = handles.elements(ind);

            %Parse out the sortable fields of the elements structure
            handles.filenums = getAnnotationVector(handles.elements,'filenum');
            handles.drugstatus = getAnnotationVector(handles.elements,'drugstatus');
            handles.directstatus = getAnnotationVector(handles.elements,'directstatus');
            
            %Set the message to show annotation replaced
            set(handles.listbox_annotation,'String',handles.annotitles)
            set(handles.text_message,'String',[file ' is added to the annotation list.']);
        end
    %This is the code for loading the first/only folder    
    else
        handles.annotationFilename = [path,file];
        handles.annotitles = file;
        annotation = aaLoadHashtable(handles.annotationFilename);

        %Copy out the two main structures of the annotation file
        handles.elements = annotation.elements;
        handles.keys = annotation.keys;

        %Sort the keys alphabetically (which is also temporally) and
        %use the index to sort the elements vector as well
        [handles.keys ind] = sort(handles.keys);
        handles.elements = handles.elements(ind);
            
        %Parse out the sortable fields of the elements structure
        handles.filenums=getAnnotationVector(handles.elements,'filenum');
        handles.drugstatus=getAnnotationVector(handles.elements,'drugstatus');
        handles.directstatus=getAnnotationVector(handles.elements,'directstatus');
        
        %Set the message to show annotation replaced
        set(handles.listbox_annotation,'String',handles.annotitles)
        set(handles.text_message,'String',[handles.annotationFilename ' is the new annotation.']);
    end
    
    %If a data location is already set, cross-reference the two
    if ~isempty(handles.totFilelist)
        %Pare down the totFilelist set to only those files that have a record
        %in the loaded annotations; dump the results in handles.filelist
        [handles.filelist] = XrefAnnotData(handles.totFilelist, handles.keys);
        set(handles.listbox_keys,'String', handles.keys)
        for i=1:length(handles.filelist)
            names{i} = handles.filelist(i).name;
            t = getFileTime(handles.filelist(i).name);
            handles.filetimes(i) = 24*(t-floor(t)); %This provides the time as the fractional hour of the day
        end
        set(handles.listbox_filelist,'String',names)
        
        %Update the record selection criteria
        %.........the options for sorting by drug status
        contents = ['All' sort(unique(handles.drugstatus))];
        set(handles.popup_drugstatus,'String',contents);
        set(handles.popup_drugstatus,'Value',1);
        handles.filter_DrugStat = contents(1);
        
        %.........the options for sorting by directed status
        contents = ['All' sort(unique(handles.directstatus))];
        set(handles.popup_directstatus,'String',contents);
        set(handles.popup_directstatus,'Value',1);
        handles.filter_DirectStat = contents(1);
        
        %.........the upper and lower bounds of the filenum selector
        set(handles.edit_fromRec,'String',num2str(min(handles.filenums)));
        handles.filter_fromRec = min(handles.filenums);
        set(handles.edit_toRec,'String',num2str(max(handles.filenums)));
        handles.filter_toRec = max(handles.filenums);
        
        %.........the upper and lower bounds of the time selector
        handles.filter_fromTime = min(handles.filetimes);
        handles.filter_toTime = max(handles.filetimes);
        set(handles.edit_startTime,'String',frac2str(handles.filter_fromTime));
        set(handles.edit_endTime,'String',frac2str(handles.filter_toTime));

    end
end
guidata(hObject, handles);

function push_clearAnnots_Callback(hObject, eventdata, handles)

button = questdlg('Are you sure you want to clear all annotations?','Clear annotations?','Clear''em All!','Nooooo!','Nooooo!');

if (strcmp(button,'Clear''em All!'))
    handles.annotationFilename = [];
    handles.annotitles = [];
    handles.filelist = [];
    handles.dataset = [];

    %Clear the two main structures of the annotation file
    handles.elements = [];
    handles.keys = [];

    %Parse out the sortable fields of the elements structure
    handles.filenums=[];
    handles.filetimes = [];
    handles.drugstatus=[];
    handles.directstatus=[];

    %Set the message to show annotation replaced
    set(handles.listbox_annotation,'String',handles.annotitles)
    set(handles.text_message,'String','All annotations cleared.');
else
    set(handles.text_message,'String','Annotations are unchanged.');
end

guidata(hObject, handles);

function push_clearFolders_Callback(hObject, eventdata, handles)
button = questdlg('Are you sure you want to clear all folders?','Clear folder?','Clear''em All!','Nooooo!','Nooooo!');

if (strcmp(button,'Clear''em All!'))
    handles.dirname = [];

    %Clear the two main repositories of the file info
    handles.totFilelist = [];
    handles.filelist = [];
    handles.dataset = [];
    
    %Set the message to show annotation replaced
    set(handles.listbox_data,'String',handles.dirname);
    set(handles.text_message,'String','All folders cleared.');
else
    set(handles.text_message,'String','Data folders are unchanged.');
end

guidata(hObject, handles);

function push_loadSylIndex_Callback(hObject, eventdata, handles)
[file,dpath]=uigetfile('*.mat');

if isequal(file,0) || isequal(dpath,0)
    set(handles.text_message,'String','Operation canceled or selected file is invalid.');
else
    load([dpath file]);
    handles.templates=templates;
    set(handles.text_message,'String',['Syllable Index was loaded from ' file '.']);
    set(handles.push_createDataset,'BackgroundColor',[0,1,0])
end
guidata(hObject, handles);

function push_createDataset_Callback(hObject, eventdata, handles)
if isempty(handles.filelist)
    set(handles.text_message,'String','No files found to load. Check annotations and data folders.')
else
    handles.dataset = [];
    h = waitbar(0,['Slow your roll, man. I got a lotta files to load here...']);
    for i=1:length(handles.filelist)
        if ~strcmp(handles.filelist(i).name,'DUMMY')
            if strcmp(handles.filelist(i).name(end),'v')
                %The file is a wav
                [handles.dataset{i}, ~] = audioread([handles.filelist(i).dirname filesep handles.filelist(i).name]);
            elseif strcmp(handles.filelist(i).name(end),'t') || strcmp(handles.filelist(i).name(end),'m')
                %The file is a dat or stm
                [chans, ~] = getChannels([handles.filelist(i).dirname filesep handles.filelist(i).name]);
                handles.dataset{i} = chans(1,:);
            end
        end
        waitbar(i/length(handles.filelist))
    end
    close(h)
    set(handles.text_message,'String','Dataset creation now completed')
end
guidata(hObject, handles);

function check_loadWAV_Callback(hObject, eventdata, handles)

function check_loadDAT_Callback(hObject, eventdata, handles)

function check_loadSTM_Callback(hObject, eventdata, handles)
%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Record Selection Controls
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_fromRec_Callback(hObject, eventdata, handles)

if str2double(get(handles.edit_fromRec,'String')) >= min(handles.filenums) && str2double(get(handles.edit_fromRec,'String')) <= str2double(get(handles.edit_toRec,'String'))
    handles.filter_fromRec = str2double(get(handles.edit_fromRec,'String'));
else
    set(handles.edit_fromRec,'String',num2str(min(handles.filenums)));
    handles.filter_fromRec = min(handles.filenums);
end
guidata(hObject, handles);

function edit_fromRec_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_toRec_Callback(hObject, eventdata, handles)

if str2double(get(handles.edit_toRec,'String')) <= max(handles.filenums) && str2double(get(handles.edit_toRec,'String')) >= str2double(get(handles.edit_fromRec,'String'))
    handles.filter_toRec = str2double(get(handles.edit_toRec,'String'));
else
    set(handles.edit_toRec,'String',num2str(max(handles.filenums)));
    handles.filter_toRec = max(handles.filenums);
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

%Limits checking
if startT >= min(handles.filetimes) && startT <= endT
    handles.filter_fromTime = startT;
else
    set(handles.edit_startTime,'String',frac2str((min(handles.filetimes))));
    handles.filter_fromTime = min(handles.filetimes);
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

%Limits checking
if endT <= max(handles.filetimes) && startT <= endT
    handles.filter_toTime = endT;
else
    set(handles.edit_endTime,'String',frac2str((max(handles.filetimes))));
    handles.filter_toTime = max(handles.filetimes);
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

%Get the remapping index
[handles.data.drugRemapInx,breaks] = sortData(handles);

%Remap the data given the produced index
sc_audio = handles.data.audioPower(handles.data.drugRemapInx,:);

handles = showAlignedData(handles,sc_audio,handles.data.drugRemapInx,breaks);

guidata(hObject, handles);

function check_sortdirect_Callback(hObject, eventdata, handles)
%Don't allow uncheck function; must select another checkbox to unselect
if ~get(handles.check_sortdirect,'Value')
    set(handles.check_sortdirect,'Value',1);
end

%Reset the linked checkboxes
set(handles.check_sortdrugs,'Value',0);
set(handles.check_sortchrono,'Value',0);

%Get the remapping index
[handles.data.directRemapInx,breaks] = sortData(handles);

%Remap the data given the produced index
sc_audio = handles.data.audioPower(handles.data.directRemapInx,:);

handles = showAlignedData(handles,sc_audio,handles.data.directRemapInx,breaks);

guidata(hObject, handles);

function check_sortchrono_Callback(hObject, eventdata, handles)
%Don't allow uncheck function; must select another checkbox to unselect
if ~get(handles.check_sortchrono,'Value')
    set(handles.check_sortchrono,'Value',1);
end

%Reset the linked checkboxes
set(handles.check_sortdrugs,'Value',0);
set(handles.check_sortdirect,'Value',0);

%Get the remapping index
[handles.data.chronoRemapInx,breaks] = sortData(handles);

%Remap the data given the produced index
sc_audio = handles.data.audioPower(handles.data.chronoRemapInx,:);

handles = showAlignedData(handles,sc_audio,handles.data.chronoRemapInx,breaks);

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File List Displays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function listbox_keys_Callback(hObject, eventdata, handles)

function listbox_keys_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_filelist_Callback(hObject, eventdata, handles)

function listbox_filelist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_buffer_Callback(hObject, eventdata, handles)

function edit_buffer_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_edgeThresh_Callback(hObject, eventdata, handles)

function edit_edgeThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_dtwBand_Callback(hObject, eventdata, handles)

function edit_dtwBand_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function check_norm_Callback(hObject, eventdata, handles)

function push_flagMotif_Callback(hObject, eventdata, handles)
%Flag the file for some sort of problem

%If no spotlight is there, generate it
if ~isfield(handles, 'spotNum')
    handles.spotNum = 1;
    
    %Update the marker and the spotlight
    handles = updateSpotlight(handles);
end

%Flip the value of the flag (i.e., if ON now OFF and opposite)
if handles.data.flag(handles.spotNum)
    handles.data.flag(handles.spotNum) = false;
else
    handles.data.flag(handles.spotNum) = true;
end

%Generate some kind of indicator so it's clear with have been flagged
updateFlags(handles);

guidata(hObject, handles);

function check_floatDTW_Callback(hObject, eventdata, handles)

function edit_offsetLims_Callback(hObject, eventdata, handles)

function edit_offsetLims_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_sigHP_Callback(hObject, eventdata, handles)

function edit_sigHP_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_sigGain_Callback(hObject, eventdata, handles)
function edit_sigGain_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_sigLP_Callback(hObject, eventdata, handles)
function edit_sigLP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function push_exportImg_Callback(hObject, eventdata, handles)
%Show data on axes
figure

[RemapIndx,breaks] = sortData(handles);
sc_audio = handles.data.audioPower(RemapIndx,:);

%Show data on axes
subplot(2,1,1)
h1=imagesc(sc_audio);
hold on
for i = 1:length(breaks)
    line([0 size(handles.data.audioPower,2)],[breaks(i) breaks(i)],'Color','w','LineWidth',1)
end
hold off

%plot basic statistics across axes_stats.  Here, I'm showing the mean
%across all trials +/- standard deviation.
subplot(2,1,2)
plot(mat2gray(mean(sc_audio)));
hold on
plot(mean(sc_audio)+std(sc_audio),':');
plot(mean(sc_audio)-std(sc_audio),':');
hold off

function push_exportRaw_Callback(hObject, eventdata, handles)

filt_keys = handles.keys(handles.filtInd==1);
filenames = filt_keys(handles.chosenStartSeq(:,1));
contents = get(handles.popup_seqSylls,'String');
sequence = contents(get(handles.popup_seqSylls,'Value'),:);
motiftimes = handles.data.motif;
audio = handles.data.audio;

[SaveName,SavePath] = uiputfile('*.mat');
save([SavePath SaveName],'sequence','filenames','motiftimes','audio','neuro');

function push_exportTemplate_Callback(hObject, eventdata, handles)

if ~isfield(handles.data,'template')
    set(handles.text_message,'String','No song template in memory. Run alignment routine and try again.');
else
    %Prep variables to save
    template = handles.data.template;
    templatesyllBreaks = handles.data.templatesyllBreaks;
    templatemotifBreaks = handles.data.templatemotifBreaks;
    
    %Get location and save file
    [fname, pathname] = uiputfile('*.mat','Select location and file name to save template.');
    save([pathname fname],'template','templatesyllBreaks','templatemotifBreaks')
    
    %Update display
    set(handles.text_message,'String',['Current template saved to file: ' pathname fname]);
end

guidata(hObject, handles);

function push_exportProc_Callback(hObject, eventdata, handles)
%Parse data to save
filt_keys = handles.keys(handles.filtInd==1);
filenames = filt_keys(handles.chosenStartSeq(:,1));
contents = get(handles.popup_seqSylls,'String');
sequence = contents(get(handles.popup_seqSylls,'Value'),:);

%Template Data
data.template = handles.data.template;
data.templatesyllBreaks = handles.data.templatesyllBreaks;
data.templatemotifBreaks = handles.data.templatemotifBreaks;

%Stuff I want
data.p = handles.data.p;
data.q = handles.data.q;
data.audio = handles.data.audio;
data.aligned_audioCube = handles.data.aligned_audioCube;

[SaveName,SavePath] = uiputfile('*.mat');
save([SavePath SaveName],'sequence','filenames','data');%,'mAudio');

%Update display
set(handles.text_message,'String',['Current process data saved to file: ' SavePath SaveName]);

function push_exportMinProc_Callback(hObject, eventdata, handles)
%Parse data to save
filt_keys = handles.keys(handles.filtInd==1);
filenames = filt_keys(handles.chosenStartSeq(:,1));
contents = get(handles.popup_seqSylls,'String');
sequence = contents(get(handles.popup_seqSylls,'Value'),:);
data.templatesyllBreaks = handles.data.templatesyllBreaks;
data.template = handles.data.template;
data.p = handles.data.p;
data.q = handles.data.q;
data.audio = handles.data.audio;
data.audioSpecs_LN = handles.data.audioSpecs_LN;
if isfield(handles.data,'aligned_audioCube')
    data.aligned_audioCube = handles.data.aligned_audioCube;
else
    data.alignedaudioPower = handles.data.audioPower;
end

data.drugsStatus = handles.data.drugsStatus;
data.directStatus = handles.data.directStatus;

[SaveName,SavePath] = uiputfile('*.mat');
save([SavePath SaveName],'sequence','filenames','data');%,'mAudio');

%Update display
set(handles.text_message,'String',['Current process data saved to file: ' SavePath SaveName]);

function push_exportPaths_Callback(hObject, eventdata, handles)

if ~isfield(handles.data,'pathR')
    set(handles.text_message,'String','No alignments paths in memory. Run alignment routine and try again.');
    return
end

%Parse data to save
paths.pathR = handles.data.pathR;
paths.pathT = handles.data.pathT;
paths.pathMean = handles.data.pathMean;
paths.p = handles.data.p;
paths.q = handles.data.q;
paths.templatesyllBreaks = handles.data.templatesyllBreaks;

[SaveName,SavePath] = uiputfile('*.mat');
save([SavePath SaveName],'paths');

%Update display
set(handles.text_message,'String',['Current paths saved to file: ' SavePath SaveName]);

function push_exportAveSpec_Callback(hObject, eventdata, handles)
if ~isfield(handles.data,'template')
    set(handles.text_message,'String','No song template in memory. Run alignment routine and try again.');
else
    %Prep variables to save
    aveSpec = squeeze(mean(handles.data.aligned_audioCube,1));
    
    %Get location and save file
    [fname, pathname] = uiputfile('*.mat','Select location and file name to save average Spectrogram.');
    save([pathname fname],'aveSpec')
    
    %Update display
    set(handles.text_message,'String',['Current average spectrogram saved to file: ' pathname fname]);
end

guidata(hObject, handles);

function push_exportTalon_Callback(hObject, eventdata, handles)
%Parse data for later use with Talon and save to file

%Retrieve key identifying information
filt_keys = handles.keys(handles.filtInd==1);
filenames = filt_keys(handles.chosenStartSeq(:,1));

filt_filenums = handles.filenums(handles.filtInd==1);
filenums = filt_filenums(handles.chosenStartSeq(:,1));

contents = get(handles.popup_seqSylls,'String');
sequence = contents(get(handles.popup_seqSylls,'Value'),:);

%Template Data
data.templatesyllBreaks = handles.data.templatesyllBreaks;
data.template = handles.data.template;

%Alignment paths
data.AlignType = get(handles.popup_alignType,'Value');
data.p = handles.data.p;
data.q = handles.data.q;
data.pFull = handles.data.pFull;
data.snipTimes = handles.data.snipTimes;

%Audio Data
audio.raw = handles.data.audio;
if isfield(handles.data,'aligned_audioCube')
    audio.aligned_audioCube = handles.data.aligned_audioCube;
else
    audio.alignedaudioPower = handles.data.audioPower;
end

%Experiment Status
data.drugsStatus = handles.data.drugsStatus;
data.directStatus = handles.data.directStatus;
data.flag = handles.data.flag;

%Request save location
[SaveName,SavePath] = uiputfile('*.mat');
save([SavePath SaveName],'sequence','filenames','filenums','data','audio');

%Update display
set(handles.text_message,'String',['Current Talon data saved to file: ' SavePath SaveName]);
