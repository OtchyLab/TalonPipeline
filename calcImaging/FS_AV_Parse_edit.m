function FS_AV_Parse_edit(source, dest)
% FS_AV_Parse_edit will do several things:
%
%   1. Seperate audio and video information from raw .MOV recordings using
%   LNP's extractmedia script
%   2. Creates new Talon-compliant filename based on raw recording names
%   3. Create/copy/move MOV, CSV, WAV, and MAT versions of the data for
%   later use
%
% FS_AVparse should be run first:
% FS_AVparse-->FS_TemplateMatch-->FS_Plot_ROI--> BatchDff2----> FS_Image_ROI--> FS_FS_Plot_ROI
%
% Parse Data from FreedomScopes
%   Created: 2015/08/02
%   By: WALIII
%   Updated: 2018/07/09 % modify file saving for use with the talon pipeline
%   By: sadieiq694
%   Updated: 2019/01/03 % streamline and reorg for real implementation
%   By: tmo

% Setup the file locations for the script
% mother = 'C:\Users\Tim\Desktop\Matlab Code\TalonPipeline\calcImaging\TestData'; %COMPUMONSTER_BU
mother = '/Users/Tim/Documents/MATLAB/TalonPipeline/calcImaging/TestData'; %MacDaddy

%Shortcut for hardcoded data source
if ~exist('source', 'var')
    source = mother;
    birdName = 'LLY78';
else
    s = regexp(source, filesep, 'split');
    birdName = s{end};
end

%Shortcut for hardcoded data destination
if ~exist('dest', 'var')
    dest = strcat(mother, filesep, birdName);
end
if exist(dest,'dir')
    rmdir(dest,'s');
end
mkdir(dest); % big directory for the bird data; need another one

% Retrieve the list of *.mov files
tmp = dir(fullfile(source,'*.mov'));
mov_listing = {tmp(:).name};
bytesList = cell2mat({tmp(:).bytes});

% Retrieve the list of *.csv files
csvfiles = dir(fullfile(source, '*.csv'));
csvfilenames ={csvfiles(:).name};

%Constants for Bandpass Audio
HP = 100; %Hz
LP = 22050; %Hz
HP_fNorm = HP/(44150/2);
LP_fNorm = LP/(44150/2);
[BP_b, BP_a] = butter(2,[HP_fNorm, LP_fNorm]);

% Cycle through the files, extracting and saving
disp('Parsing Audio and Video files');
for i=1:length(mov_listing)
    
    % If not existent, create directory based on date of .mov file
    date = mov_listing{i}(1:10);
    specTarget = strcat(dest, filesep, date);
    if ~exist(specTarget, 'dir')
        mkdir(specTarget); 
    end

    %Current working file; get file details
    curFile = [source, filesep, mov_listing{i}];
    
    %To prevent crashes/long processing time, restrict processing to
    %shorter recordings (i.e., <900 MB)
    if bytesList(i)/1000000< 900.000
        %Extract Audio and video data...
        try
            %MacOS-only extraction code
            [a_ts, a, v_ts, v] = extractmedia(curFile);
            
        catch
            %This does not appear to work like the extractmedia routine
            %above... needs troubleshooting
            v1 = VideoReader(curFile);
            k = 1;
            while hasFrame(v1)
                v{k} = readFrame(v1);
                k = k+1;
            end
            v = v';
            a = 0;
            a_ts = 0;
            v_ts = 0;
            
        end
        
    else
        
        disp('Current file too large for batch processing... check the function limits');
        %LargeDir = strcat(path,'/','LargeFiles');
        %movefile(FILE, error_dir)
        continue;
    end

    %%%%%%%%%%%%%%%%%%%
    % Construct video block
    %%%%%%%%%%%%%%%%%%%
    
    % Extract video format
    video = [];
    [video.width, video.height, video.channels] = size(v{1});
    video.times = v_ts; % 0.1703*(day-1);
    video.nrFramesTotal = size(v,1);
    video.FrameRate = 1/mean(diff(v_ts));
    
    %Stack the cells as frames from each RGB channel
    est = [];
    for ii = 1:size(v,1)
        est(:,:,:,ii) = v{ii};
    end
    video.frames = est;
    
    %Noise calculations -- don't know that I really buy these numbers...
    noise = squeeze(est(:,:,3,:)); % blue channel (no imaging signals here)
    noise = (((squeeze(mean(mean(noise(:,1:80,:),1))))));%+(squeeze(mean(mean(noise(1:10,:,:),2)))))/2);
    video.gain = noise;
    
    %%%%%%%%%%%%%%%%%%%
    % Construct audio block
    %%%%%%%%%%%%%%%%%%%
    
    % Change audio data to proper sampling frequency
    a_filt = filtfilt(BP_b, BP_a, double(a));
    [P,Q] = rat(44150/48000);
    rs_audio = resample(a_filt, P, Q); %resample with FIR (better than simple interpolation)
    
    %Audio output structure
    audio = [];
    audio.nrChannels = 1;
    audio.bits = 16;
    audio.nrFrames = length(rs_audio);
    audio.data = rs_audio;
    audio.rate = 44150;
    audio.TotalDurration = audio.nrFrames/audio.rate;
    
    %%%%%%%%%%%%%%%%%%%
    %Renaming/sorting
    %%%%%%%%%%%%%%%%%%%
    
    % Create file timestamp serial number
    [Y, M, D, H, MN, S] = getDateVec(mov_listing{i}); %Parse filename
    timeStamp = datenum([str2num(Y), str2num(M), str2num(D), str2num(H), str2num(MN), str2num(S)]);
    refTime = datenum([1903 12 31 20 00 00]);
    serialNum = num2str((timeStamp-refTime)*(24*60*60)+1);
    
    % Assemble filename base from components
    name = formFileName(birdName, serialNum, Y, M, D, H, MN, S);
    longname = [specTarget, filesep, char(name)];
    
    % Write/move files to new target
    copyfile(mov_listing{i}, [longname, '.mov']);
    copyfile(csvfilenames{i}, [longname, '.csv']);
    audiowrite([longname, '.wav'], audio.data, audio.rate);
    save([longname, '.mat'], 'audio', 'video', '-v7.3');
    
end

