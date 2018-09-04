function [start_Syl, end_Syl] = old_seqSyl(audio, min_delai_threshold, min_t_threshold,ratio)
% old_seqSyl  separate an audio file into subparts of interest using the
% power spectrum, one amplitude threshold and two temporal thresholds. 
% 
% AUDIO                         is the audio file
% MIN_DELAI_THRESHOLD           is the threshold (in number of points) for
% the delai between two consecutive interesting parts of the audio file
% (two consecutive parts above the amplitude threshold). If the delai
% between these two regions is below this value, we consider the two
% regions as one. It is set as 50ms in my routine process. 
% MIN_T_THRESHOLD               is the temporal threshold to eliminate the
% smaller intereesting part of the file. It is set as
% (length_minimal_song_file-(2*0.5*44150))/2 in my routine process. This is
% half of the length of the minimal song file without the silent parts.
% RATIO                         is the amplitude threshold above which the
% audio file is interesting. The threshold is set by RATIO time the sd of
% background noise. It can be consider as a signal/noise ratio.

%set the amplitude threshold and apply it to the audio file
amplitude_threshold=ratio*eval_bg_noise(audio,50,20000);
X=vertcat(zeros(49,1), ((audio(50:end).^2)>amplitude_threshold));X(1)=0;X(end)=0;

diff1=find(diff(X)==1); %start of each period of interest
diff_1=find(diff(X)==-1); %end of each period of interest
delai=vertcat(diff1,0) - vertcat(0,diff_1); % delai between two consecutive regions

% applying min_delai threshold to suppress two consecutive regions which
% are too close in time.
vector_del=find(delai<min_delai_threshold);vector_del(end)=[];
if vector_del(1)==1
    vector_del(1)=[];
end
diff1(vector_del)=[];
diff_1(vector_del-1)=[];

% applying the second time threshold to consider only the interesting
% region longer than this threshold
if isempty(diff1)==0
start_Syl=diff1((diff_1-diff1)>min_t_threshold);
end_Syl=diff_1((diff_1-diff1)>min_t_threshold);
end

end
