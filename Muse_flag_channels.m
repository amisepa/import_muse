
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_bad = .33;  % max fraction of windows to be tolerated before flagging channel
win_size = 5;   % window size (in s) CLASSIFIERS TRAINED FOR 5-s WINDOWS
usegpu = 0;     % use GPU computing
vis = 1;        % visualize flagged channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load classifiers (one for frontal and one for posterior because they have 
% pretty different signals)
mdlPath = 'add_path_to_classifiers'; %%%%%%%%% ADD PATH %%%%%%%%
load(fullfile(mdlPath, 'trainedModelFront\trainedModelFrontLogisticReg.mat'));  % default = trainedModelFrontLogisticReg.mat
load(fullfile(mdlPath, 'trainedModelPost\trainedModelPostTree.mat'));           % default = trainedModelPostTree.mat

%%%%%%%%% Load a MUSE file %%%%%%%%

% Filter (DO NOT EDIT FOR CLASSIFIERS - NOT TESTED WITH OTHER CUTOFF FREQS)
EEG = pop_eegfiltnew(EEG,'locutoff',1);    
EEG = pop_eegfiltnew(EEG,'hicutoff',50);   

% Scan channels to flag them as bad
[badChan, badChanLabels] = scan_channels(EEG, trainedModelFront, trainedModelPost, max_bad, win_size, usegpu, vis);
