function audioData = changeFrameRate(raw, rate_orig, rate_new)
% fixing the frequency of audio recordings! :)

%Size and content of the raw data vector
x = 1:numel(raw);
v = raw';

%Number of datapoints at the new sampling rate
n = round((numel(raw)/rate_orig) * rate_new);

%Lookup times
xq = linspace(1, numel(v), n);

%resample via simple interpolation
audioData = interp1(x, v, xq);