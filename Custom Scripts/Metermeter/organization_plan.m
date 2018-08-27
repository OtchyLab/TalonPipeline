% Three types of analysis: Variance Decomposition, Spectral Analysis,
% Temporal Analysis
% Create different functions for each part of the process (which all of
% them go through), call each of these functions in the appropriate button
% callbacks
% 1) Calculations
% 2) Plotting 
% 3) Saving
% Ideally each callback would call one function; right now some of the work
% is done in the callbacks themselves. Change this once everything's
% working smoothly

% Variance Decomposition
% 1) Calculations - can be found in the push_varDecomp_Callback
% 2) Plotting - plotVarComp -> called by push_plotVarianceDecomp_Callback
% 3) Saving - can be found in the push_varianceDecompSave_Callback
%       Not sure this is saving exactly what is wanted/required yet

% Spectral Analysis
%   Breaks down to syllable similarity and syllable recovery
%       Syllable similarity:
%       - Calculations: calculations found in push_syllSim_Callback
%       - Plotting: plotSimAnalysis with 1 as the caller (called in
%       push_plotSyllSim_Callback)
%       - Saving: can be found in push_spectralAnalysisSave_Callback; saves
%       both similarity and recovery if they exist
%       Syllable Recovery: 
%       - Calculations:
%       - Plotting: plotSimAnalysis with 0 as the caller (called by
%       push_plotSyllRecovery_Callback)
%       - Saving: can be found in push_spectralAnalysisSave_Callback; saves
%       both similarity and recovery if they exist

% Temporal Analysis
% 1) Calculations: 
% 2) MAKE ONE PLOT FUNCTION
function handles = plotTempAnalysis(handles)
% get the value in the dropdown menu
plotType = get(handles.popup_

% Tackle temporal analysis later! Scary and you're tired. You have to
% delete a lot of buttons and understand what calculations are done that
% get saved. Also get the plot to change only when you press the plot
% button and the type of plot to depend on the entry in the dropdown menu 


% Also -- find a better way to select intervals? And do intervals for
% syllable recovery? I forget how they are selected for that. Anyway, the
% current method of interval selection is rather effective as it puts no
% limits on # of intervals, which most other options would. Ask how
% important it is to be able to select multiple intervals/what the usual
% number of intervals that is selected is before making a decision


