function min_sd_noise=eval_bg_noise(audio, begin_noise, end_noise)
% eval_bg_noise     permit to approximate the standart deviation of the
% background noise. The function divide an part of a wav file (specify by
% the input) in 10 egal parts and estimate the sd on each one. The output
% is the minimal value of the sd obtained. It procced like this to avoid
% biais due to high noise.

% AUDIO             is the audio file imported by WAVREAD
% BEGIN_NOISE       is the beginning of the part of the audio file we want
% to treat (in number of points) 
% END_NOISE         is the end of the part of the audio file we want to
% treat (in number of points) 

% MIN_SD_NOISE      is the smaller sd obtain on the 10 subparts. 

all_sd_noise=zeros(10,1);
index=linspace(begin_noise,end_noise,10);
parfor i=1:9
    all_sd_noise(i)=sqrt(var(audio(floor(index(i)):floor(index(i+1))).^2));
end
min_sd_noise=min(all_sd_noise);
end