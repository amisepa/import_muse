%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Usage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flag bad Muse EEG channels using the trained classifiers shipped with the
% import_muse plugin (classifier_front.mat + classifier_post.mat, loaded
% automatically by scan_channels). EEG must be a raw (unprocessed) Muse EEG
% structure imported with import_muse (4 channels: TP9, AF7, AF8, TP10).
%
%   EEG = import_muse(file_path);                          % import raw Muse EEG
%   EEG = pop_eegfiltnew(EEG,'locutoff',1,'hicutoff',50);  % filter 1-50 Hz
%   [badChan, badChanLabels] = scan_channels(EEG, .33, 1); % flag channels
%
% Or in one step during import:
%   EEG = import_muse(file_path, 'detectBadChan');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_bad = .33;  % max fraction of bad 5-s windows tolerated before flagging the channel
win_size = 5;   % window size (in s) CLASSIFIERS TRAINED FOR 5-s WINDOWS
usegpu = 0;     % use GPU computing (unused)
vis = 1;        % visualize flagged channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter (DO NOT EDIT FOR CLASSIFIERS - NOT TESTED WITH OTHER CUTOFF FREQS)
% NOTE: scan_channels() also filters on the fly; filter here only if calling
% scan_channels yourself on unfiltered data.
EEG = pop_eegfiltnew(EEG,'locutoff',1);
EEG = pop_eegfiltnew(EEG,'hicutoff',50);

% Scan channels to flag them as bad
[badChan, badChanLabels] = scan_channels(EEG, max_bad, vis);
