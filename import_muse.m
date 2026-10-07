%% import_muse() - Import Muse .csv data recorded with the Mind Monitor,
% Muse Direct, or Muse Lab (OpenMuse) Apps, for all Muse headsets
% (Muse 1, Muse 2, Muse S, and Muse S Athena incl. its fNIRS/optics data).
% EEG data are always imported by default. Optional inputs include ACC, GYR,
% AUX, PPG, and OPTICS channels.
%
% Requirements: MATLAB and EEGLAB.
%
% Input:
%   file_path   - full path and name of the Muse file (e.g. 'file_path\file_name.csv')
%
% Optional inputs:
%   'acc'            - Import accelerometer data (default OFF = 0)
%   'gyr'            - Import gyroscope data (default OFF = 0)
%   'ppg'            - Import photoplethysmogram data (default OFF = 0; Muse Direct recordings only)
%   'aux'            - Import auxiliary channel data (default OFF = 0; MindMonitor recordings only)
%   'optics'         - Import optical fNIRS data (default OFF = 0; Muse S Athena only)
%   'detectBadChan'  - Detect and remove bad channels using the trained classifiers
%
% Outputs:
%   EEG     - Data in EEGLAB structure format containing signal from each
%             selected channel (time-synchronized). If the file contains no
%             EEG channels (e.g. an optics-only recording), the structure
%             contains the optical channels instead.
%   com     - history string for command line.
%
% Usage:
%   EEG = import_muse;                                          % select file and parameters in GUI mode
%   EEG = import_muse(file_path);                               % import EEG data using command line
%   EEG = import_muse(file_path, 'ppg');                        % import EEG and PPG signal
%   EEG = import_muse(file_path, 'acc', 'gyr', 'aux', 'ppg');   % import everything
%   EEG = import_muse(file_path, 'optics');                     % import EEG and fNIRS optical signal
%   EEG = import_muse(file_path, 'detectBadChan');              % import EEG and flag/remove bad channels
%
% Important notes:
%   ACC and GYR amplitude is modified to better fit the EEG data, for plotting
%   purposes; see command line for correction. Athena optics (fNIRS) channels
%   are imported raw (no amplitude modification) at their 64 Hz recording rate,
%   or resampled to the EEG rate when requested together with EEG signals.
%
% Reference: Cannard, C., Wahbeh, H., & Delorme, A. (2021). Validating the
% wearable MUSE headset for EEG spectral analysis and Frontal Alpha Asymmetry.
% 2021 IEEE International Conference on Bioinformatics and Biomedicine (BIBM),
% 3603-3610. https://doi.org/10.1109/BIBM52615.2021.9669778
%
% Author: Cedric Cannard, CerCo, CNRS
%
% Copyright (C) Cedric Cannard, 2020-2026
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version. This program is distributed in
% the hope that it will be useful, but WITHOUT ANY WARRANTY, without even
% the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG, com] = import_muse(file_path, varargin)

%% File path and name using pop-up window
if nargin < 1
    [fileName, filePath] = uigetfile2({'*.csv'}, 'Select Muse .csv file - import_muse()');
    if isempty(fileName) || isequal(fileName, 0)
        EEG = []; com = '';
        disp('No file selected; import cancelled.');
        return;
    end
    file_path = fullfile(filePath, fileName);
end

%% IMPORT DATA
disp('Importing EEG data...');
data = readtable(file_path, 'Delimiter', ',', 'VariableNamingRule', 'preserve');

% Read the header line separately (robust across MATLAB versions, unlike importdata)
fID = fopen(file_path, 'r');
header = fgetl(fID);
fclose(fID);
if ~isempty(header)
    % strip byte-order marks (UTF-8 BOM = U+FEFF char or EF BB BF bytes mis-decoded)
    while ~isempty(header) && (double(header(1)) == 65279 || double(header(1)) < 9 || (double(header(1)) > 13 && double(header(1)) < 32))
        header(1) = [];
    end
end
if isempty(header)
    error('Could not read the header of %s (empty or corrupted file?)', file_path);
end
varNames = lower(strsplit(header, ','));   % cellstr, same order as table columns

%% FORMAT DETECTION
% Mind Monitor files:        'TimeStamp',Delta_TP9,...,RAW_TP9,...    (datetime/mm:ss.S text)
% Muse Direct files:         'timestamps',eeg_1,...                   (POSIX seconds)
% Athena (Muse Direct app):  'Timestamp','PacketType','Data'          (POSIX microseconds; one packet per row)
% Athena (Muse Lab app):     'timestamp','osc_address','osc_type','osc_data'  (POSIX milliseconds; OSC log)
% Athena (EEG export):       'ts',TP9,AF7,AF8,TP10                    (POSIX seconds)
% Athena (optics export):    'ts',ch1,...,ch16                        (POSIX seconds; 64 Hz)
fileFormat = 'classic';
if numel(varNames) >= 2 && strcmpi(varNames{2}, 'packettype')
    fileFormat = 'athena_packets';
    rec_type = 'muse_direct_athena';
    disp('Recording App detected: Muse Direct (Muse S Athena packets)');
elseif numel(varNames) >= 2 && strcmpi(varNames{2}, 'osc_address')
    fileFormat = 'athena_osc';
    rec_type = 'muselab_athena';
    disp('Recording App detected: Muse Lab / OpenMuse OSC log (Muse S Athena)');
elseif strcmpi(varNames{1}, 'ts') && any(strcmpi(varNames, 'tp9'))
    fileFormat = 'athena_eeg';
    rec_type = 'muse_direct_athena';
    disp('Recording App detected: Muse Direct (Muse S Athena EEG export)');
elseif strcmpi(varNames{1}, 'ts') && any(startsWith(varNames, 'ch'))
    fileFormat = 'athena_optics';
    rec_type = 'muse_direct_athena';
    disp('Recording App detected: Muse Direct (Muse S Athena optics/fNIRS export)');
elseif ~contains(varNames{1}, 'time')
    error('First column of the file must contain the timestamps (got ''%s''); this does not look like a Muse recording file', varNames{1});
end

if strcmpi(fileFormat, 'classic')
    %% TIMESTAMPS (classic files: Mind Monitor and Muse Direct)
    Time = data{:,1};
    if isdatetime(Time)                          % MindMonitor file (full datetime)
        rec_type = 'muse_monitor';
        disp('Recording App detected: MindMonitor');
        Time = datetime(Time, 'Format', 'HH:mm:ss.SSS');
    elseif iscell(Time)                          % MindMonitor file (mm:ss.S strings)
        rec_type = 'muse_monitor';
        disp('Recording App detected: MindMonitor');
        % legacy clock strings (mm:ss.S) roll over at 60 min: make monotonically increasing first
        Time = datetime(Time, 'InputFormat', 'mm:ss.S', 'Format', 'HH:mm:ss.SSS');
        Time = Time - Time(1);                   % durations from recording start
        wraps = cumsum([0; diff(seconds(Time)) < -1800]);  % backward jump > 30 min = hour rollover
        Time = datetime('today') + seconds(seconds(Time) + wraps * 3600);  % unwrap rollovers + back to a datetime grid
        Time = Time(:);
    elseif isnumeric(Time)                       % Muse Direct file (POSIX seconds)
        rec_type = 'muse_direct';
        disp('Recording App detected: Muse Direct');
        Time = datetime(Time, 'ConvertFrom', 'posixtime', 'Format', 'HH:mm:ss.SSS');
    else
        error('Timestamps format not recognized! Please submit an issue on Github and attach your raw file.');
    end

    %% EEG DATA (classic files)
    % MindMonitor raw EEG columns are named RAW_TP9, RAW_AF7, RAW_AF8, RAW_TP10.
    % Muse Direct raw EEG columns are named eeg_1 to eeg_6 (first 4 = raw EEG).
    ind_eeg = find(startsWith(varNames, 'raw_'));   % header indices (col 1 = time)
    if isempty(ind_eeg)
        ind_eeg = find(startsWith(varNames, 'eeg_'));
        if numel(ind_eeg) > 4
            ind_eeg = ind_eeg(1:4);
            warning('More than 4 raw EEG columns found; importing the first 4 (eeg_1 to eeg_4).');
        end
    end
    if isempty(ind_eeg)
        error('No raw EEG columns found in the file (RAW_xx or eeg_N missing)');
    end
    ind_eeg = ind_eeg - 1;   % convert to table column indices (time column removed)
    data(:,1) = [];
    eegData = table2array(data(:,ind_eeg));

    % Remove empty rows (rows without EEG samples, e.g. band power, markers, metadata)
    nans = isnan(eegData(:,1));
    eegData(nans,:) = [];
    Time(nans,:) = [];
    nEEG = size(eegData,1);
    if nEEG < 10
        error('No raw EEG samples found in the file (only %d EEG row(s))', nEEG);
    end

    %% SAMPLE RATE (mode of EEG sample counts per second of recording)
    eeg_sRate = local_sRateFromTime(Time, nEEG);

    % Summary/band-power files contain about 1 row per second: no raw EEG
    if eeg_sRate < 50
        error(['This file contains no raw EEG samples (about ' num2str(eeg_sRate, '%g') ' rows/s). ', ...
               'It only holds summary/band-power rows (or the recording failed). Raw EEG export must be enabled in the recording app.']);
    end

    % If sample rate is too far (50 Hz) from the manufacturer default (256 Hz), cancel all actions (i.e., bad file)
    if abs(256-eeg_sRate) > 50
        % diagnose the common cause: unreliable timestamps (overlapping/duplicate rows)
        tRaw = datenum(Time);
        nBack = sum(diff(tRaw) < 0);
        if nBack > 0
            error(['Sample rate = %g Hz --> far from manufacturer''s default (256 Hz). ', ...
                   'The timestamps are not monotonic (%d backwards steps: duplicate/overlapping rows). ', ...
                   'This recording is corrupted (Mind Monitor wrote overlapping second-blocks); re-record.'], eeg_sRate, nBack);
        end
        error('Sample rate = %g Hz --> far from manufacturer''s default (256 Hz)', eeg_sRate);
    end
    eegLabels = {'TP9' 'AF7' 'AF8' 'TP10'};
end %% end classic-file parsing

%% Optional inputs (other signals)

params = struct('acc', 0, 'gyr', 0, 'ppg', 0, 'aux', 0, 'optics', 0);

%GUI
if nargin < 1
    uilist = {
        {'Style' 'text' 'string' 'Which other signals do you wish to import?'} ...
        {'style' 'checkbox' 'string' 'Accelerometer (ACC)' 'tag' 'acc' 'value' 0 'enable' 'on' } ...
        {'style' 'checkbox' 'string' 'Gyroscope (GYR)' 'tag' 'gyr' 'value' 0 'enable' 'on' } ...
        {'style' 'checkbox' 'string' 'Import Photoplethysmogram (PPG; Muse 2 and S recorded with Muse Direct only)' 'tag' 'ppg' 'value' 0 'enable' 'on' } ...
        {'style' 'checkbox' 'string' 'Import Auxiliary (AUX; Muse 1 recorded with MindMonitor only)' 'tag' 'aux' 'value' 0 'enable' 'on' } ...
        {'style' 'checkbox' 'string' 'Import optical fNIRS data (OPTICS; Muse S Athena only)' 'tag' 'optics' 'value' 0 'enable' 'on' } ...
        {} ...
        };
    uigeom = { 1 1 1 1 1 1 };
    opt = inputgui(uigeom, uilist, 'pophelp(''import_muse'')', ['Muse data recorded with ' rec_type]);
    if isempty(opt)
        EEG = []; com = '';   % user cancelled
        disp('Import cancelled.');
        return;
    end
    params.acc = opt{1};
    params.gyr = opt{2};
    params.ppg = opt{3};
    params.aux = opt{4};
    params.optics = opt{5};
else
    opt = varargin;
    for iOpt = 1:length(opt)
        switch lower(opt{iOpt})
            case 'acc', params.acc = 1;
            case 'gyr', params.gyr = 1;
            case 'ppg', params.ppg = 1;
            case 'aux', params.aux = 1;
            case 'optics', params.optics = 1;
            case 'detectbadchan'  % handled at the end of the import
            otherwise
                warning('Unknown option ''%s'' ignored (use acc, gyr, ppg, aux, optics, detectBadChan)', string(opt{iOpt}));
        end
    end
end

%% ATHENA FILES (packet, OSC, EEG-only, optics-only)

if strcmpi(fileFormat, 'athena_eeg')
    % Wide EEG-only export: ts,TP9,AF7,AF8,TP10 (POSIX seconds)
    eegData = table2array(data(:, 2:5));
    Time = datetime(data{:,1}, 'ConvertFrom', 'posixtime', 'Format', 'HH:mm:ss.SSS');
    eegData = double(eegData);
    eeg_sRate = local_sRateFromTime(Time, size(eegData,1));
    eegLabels = {'TP9' 'AF7' 'AF8' 'TP10'};
    nans = false(size(eegData,1), 1);   % no empty rows in this format

elseif strcmpi(fileFormat, 'athena_optics')
    % Wide optics-only export: ts,ch1,...,ch16 (64 Hz optical data)
    opticData = table2array(data(:, 2:end));
    Time = datetime(data{:,1}, 'ConvertFrom', 'posixtime', 'Format', 'HH:mm:ss.SSS');
    opticData = double(opticData);
    nRows = size(opticData,1);
    if nRows < 10, error('No optical samples found in the file'); end
    opt_sRate = local_sRateFromTime(Time, nRows);
    if opt_sRate < 20 || opt_sRate > 200
        error('Sample rate = %g Hz, not a valid optic rate (~64 Hz)', opt_sRate);
    end
    opticLabels = local_opticLabels(size(opticData,2));
    disp(['Optical (fNIRS) sample rate detected: ' num2str(opt_sRate) ' Hz, ' num2str(size(opticData,2)) ' channels']);
    % Return the optical data as the main structure (this file has no EEG)
    EEG = eeg_emptyset;
    EEG.chanlocs = struct('labels', {opticLabels{:}});
    EEG.data = opticData';
    EEG.srate = opt_sRate;
    EEG.pnts = size(EEG.data,2);
    EEG.nbchan = size(EEG.data,1);
    EEG.xmin = 0;
    EEG.trials = 1;
    EEG.setname = 'Optical fNIRS data (Muse S Athena)';
    EEG = eeg_checkset(EEG);
    disp('Athena optical (fNIRS) data were imported into EEGLAB.');
    com = sprintf('EEG = import_muse(''%s'', ''optics'');', file_path);
    return;

elseif any(strcmpi(fileFormat, {'athena_packets', 'athena_osc'}))
    % One packet/OSC message per row; find EEG and optic rows
    if strcmpi(fileFormat, 'athena_packets')
        pt = lower(string(data.("PacketType")));
        payload = data.("Data");
        isEEG = pt == 'eeg';
        isOpt = pt == 'optics';
        isACC = pt == 'accelerometer';
        isGYR = pt == 'gyro';
        tSec = double(data.("Timestamp")) / 1e6;    % POSIX microseconds -> seconds
    else
        addr = lower(string(data.("osc_address")));
        payload = data.("osc_data");
        isEEG = addr == '/henderson_lab/eeg';
        isOpt = addr == '/henderson_lab/optics';
        isACC = addr == '/henderson_lab/acc';
        isGYR = addr == '/henderson_lab/gyro';
        tSec = double(data.("timestamp")) / 1e3;    % POSIX milliseconds -> seconds
    end

    % EEG packets (4 raw EEG values per packet, 256 Hz; AUX_L/AUX_R on newer presets may add more)
    eegRows = find(isEEG);
    if isempty(eegRows) && ~params.optics
        error('No EEG packets found in this file; use the ''optics'' input to import the optical fNIRS data');
    end
    eegData = [];
    Time = [];
    if ~isempty(eegRows)
        [eegData, Time] = local_parseRows(payload, tSec, eegRows);
        nEEG = size(eegData,1);
        if nEEG < 10
            error('No raw EEG samples found in the file (only %d EEG row(s))', nEEG);
        end
        eeg_sRate = local_sRateFromTime(Time, nEEG);
        if eeg_sRate < 50 || abs(256-eeg_sRate) > 50
            error('EEG sample rate = %g Hz, far from manufacturer''s default (256 Hz)', eeg_sRate);
        end
        eegLabels = {'TP9' 'AF7' 'AF8' 'TP10'};
        if size(eegData,2) > numel(eegLabels)
            eegLabels{end+1:size(eegData,2)} = arrayfun(@(k) sprintf('AUX%d', k), 1:size(eegData,2)-numel(eegLabels), 'uni', 0); %#ok<AGROW>
        end
    elseif params.optics && isempty(find(isOpt, 1))
        error('No EEG and no OPTICS packets found in this file (nothing to import)');
    end

    % Optical packets (Athena fNIRS; 64 Hz; 4, 8, or 16 values per packet)
    optRows = find(isOpt);
    opticData = []; optTimenum = [];
    if ~isempty(optRows)
        [opticData, optSec] = local_parseRows(payload, tSec, optRows);
        if size(opticData,1) >= 10
            opt_sRate = local_sRateFromTime(datetime(optSec, 'ConvertFrom', 'posixtime'), size(opticData,1));
            opticLabels = local_opticLabels(size(opticData,2));
            disp(['Optical (fNIRS) data: ' num2str(size(opticData,1)) ' samples @ ' num2str(opt_sRate) ' Hz, ' num2str(size(opticData,2)) ' channels']);
            optTimenum = optSec;
        else
            opticData = [];
            disp('Fewer than 10 optical samples found; skipping OPTICS.');
        end
    end

    % Accelerometer packets
    accData = []; accTimenum = [];
    if ~isempty(find(isACC, 1))
        [accData, accSec] = local_parseRows(payload, tSec, find(isACC));
        accTimenum = accSec;
    end

    % Gyroscope packets
    gyroData = []; gyroTimenum = [];
    if ~isempty(find(isGYR, 1))
        [gyroData, gyroSec] = local_parseRows(payload, tSec, find(isGYR));
        gyroTimenum = gyroSec;
    end

    if isempty(eegData) && isempty(opticData)
        error('No EEG and no OPTICS samples found in this file (nothing importable)');
    end
    if isempty(eegData) && ~isempty(opticData)
        % Optics-only packet file: return optical data as the main structure
        EEG = eeg_emptyset;
        EEG.chanlocs = struct('labels', {opticLabels{:}});
        EEG.data = opticData';
        EEG.srate = opt_sRate;
        EEG.pnts = size(EEG.data,2);
        EEG.nbchan = size(EEG.data,1);
        EEG.xmin = 0;
        EEG.trials = 1;
        EEG.setname = 'Optical fNIRS data (Muse S Athena)';
        EEG = eeg_checkset(EEG);
        disp('Athena optical (fNIRS) data were imported into EEGLAB.');
        optFlags = varargin(~strcmpi(varargin, 'eeg'));
        if isempty(optFlags)
            com = sprintf('EEG = import_muse(''%s'', ''optics'');', file_path);
        else
            com = sprintf('EEG = import_muse(''%s'', %s);', file_path, vararg2str(optFlags));
        end
        return;
    end
    nans = false(size(eegData,1), 1);   % no empty rows in these formats
else
end

%% EEG

if strcmpi(fileFormat, 'athena_optics')
    % optics-only file already returned above
else
    EEG = eeg_emptyset;
    EEG.chanlocs = struct('labels', {eegLabels{:}});
    EEG.data = double(eegData');
    EEG.srate = eeg_sRate;
    EEG.pnts   = size(EEG.data,2);
    EEG.nbchan = size(EEG.data,1);
    EEG.xmin = 0;
    EEG.trials = 1;
    EEG.setname = ['EEG data (' rec_type ')'];
    EEG = eeg_checkset(EEG);

    % helper time grid (datenum) for resampling the ancillary signals
    if exist('Time', 'var')
        eegTimenum = datenum(Time);
    end
    if ~exist('eegTimenum', 'var') || isempty(eegTimenum)
        eegTimenum = (0:EEG.pnts-1)'/EEG.srate;
    end
end

%% OPTICS (Athena fNIRS): resample to the EEG time grid and append when requested

if params.optics && (exist('opticData', 'var') && ~isempty(opticData)) && exist('EEG','var')
    disp('Importing OPTICS (fNIRS) data.');
    if size(opticData,2) ~= numel(opticLabels)
        opticLabels = local_opticLabels(size(opticData,2));
    end
    if size(opticData,1) == EEG.pnts
        optResampled = opticData;
    else
        disp('Resampling OPTICS data to match EEG sampling rate...');
        % duplicate timestamps (ms-truncated) break interp1: keep the first row per unique time
        [uT, ia] = unique(optTimenum(:));
        optResampled = interp1(uT, opticData(ia,:), eegTimenum(:), 'nearest', 'extrap');
    end
    for iO = 1:size(optResampled,2)
        nChans = size(EEG.data,1);
        EEG.data(nChans+1,:) = optResampled(:,iO)';
        EEG.chanlocs(nChans+1).labels = opticLabels{iO};
    end
    EEG.nbchan = size(EEG.data,1);
    EEG = eeg_checkset(EEG);
elseif params.optics && (~exist('opticData', 'var') || isempty(opticData))
    warning('No optical (fNIRS) data found in this file; skipping OPTICS.');
end

%% ACC

if params.acc
    disp('Importing ACC data.');
    if any(strcmpi(fileFormat, {'athena_packets','athena_osc'}))
        if isempty(accData), warning('No accelerometer packets found in this file.'); else
            if size(accData,1) == EEG.pnts
                accRes = accData;
            else
                % duplicate timestamps (ms-truncated) break interp1: unique rows only
                [uT, ia] = unique(accTimenum(:));
                accRes = interp1(uT, accData(ia,:), eegTimenum(:), 'nearest', 'extrap');
            end
            % Transform ACC data to match EEG amplitude (for plotting purposes)
            amp_acc = std(accRes,1);   % amplitude (gravity DC offset excluded)
            amp_eeg = mean(std(EEG.data,[],2));
            d = (amp_eeg ./ amp_acc) / 3;
            accRes = accRes .* d;
            disp(['ACC amplitudes were multiplied by ' num2str(d(1)) ' (X), ' num2str(d(2)) ' (Y), ' num2str(d(3)) ' (Z) to match EEG data scale.']);
            nChans = size(EEG.data,1);
            EEG.data(nChans+1:nChans+3,:) = accRes';
            EEG.chanlocs(nChans+1).labels = 'ACC_X';
            EEG.chanlocs(nChans+2).labels = 'ACC_Y';
            EEG.chanlocs(nChans+3).labels = 'ACC_Z';
            EEG.nbchan = EEG.nbchan+3;
            EEG = eeg_checkset(EEG);
        end
    else
        ind_acc = find(startsWith(varNames, 'acc_') | startsWith(varNames, 'accelerometer_'));
        ind_acc = ind_acc - 1;
        accData = table2array(data(:,ind_acc));

        if strcmp(rec_type, 'muse_monitor')     %same sampling rate as EEG
            accData(nans,:) = [];
        else                                    %muse_direct: different sampling rate
            disp('Filling empty values of ACC data to match EEG sampling rate...');
            vals = find(~isnan(accData(:,1)));
            for i = 1:size(vals,1)
                Vals = accData(vals(i),:);
                if i == 1
                    for j = 1:vals(i)
                        accData(j,:) = Vals;
                    end
                else
                    for j = vals(i-1)+1:vals(i)
                        accData(j,:) = Vals;
                    end
                end
                if i == size(vals,1)
                    for j = vals(i):size(accData,1)
                        accData(j,:) = Vals;
                    end
                end
            end
            accData(nans,:) = [];    %Remove NaNs to match EEG data length
        end

        % Transform ACC data to match EEG amplitude
        amp_acc = std(accData,1);   % amplitude (gravity DC offset excluded)
        amp_eeg = mean(std(EEG.data,[],2));
        for i = 1:length(amp_acc)
            d(i) = (amp_eeg / amp_acc(i)) / 3 ; %#ok<AGROW>
            accData(:,i) = bsxfun(@times, accData(:,i), d(i));
        end
        disp(['ACC_X amplitude was multiplied by ' num2str(d(1)) ' to match EEG data scale.']);
        disp(['ACC_Y amplitude was multiplied by ' num2str(d(2)) ' to match EEG data scale.']);
        disp(['ACC_Z amplitude was multiplied by ' num2str(d(3)) ' to match EEG data scale.']);

        %Add to EEG structure
        nChans = size(EEG.data,1);
        EEG.data(nChans+1:nChans+3,:) = accData';
        EEG.chanlocs(nChans+1).labels = 'ACC_X';
        EEG.chanlocs(nChans+2).labels = 'ACC_Y';
        EEG.chanlocs(nChans+3).labels = 'ACC_Z';
        EEG.nbchan = EEG.nbchan+3;
        EEG = eeg_checkset(EEG);
    end
end

%% GYRO

if params.gyr
    disp('Importing GYRO data.');
    if any(strcmpi(fileFormat, {'athena_packets','athena_osc'}))
        if isempty(gyroData), warning('No gyroscope packets found in this file.'); else
            if size(gyroData,1) == EEG.pnts
                gyrRes = gyroData;
            else
                % duplicate timestamps (ms-truncated) break interp1: unique rows only
                [uT, ia] = unique(gyroTimenum(:));
                gyrRes = interp1(uT, gyroData(ia,:), eegTimenum(:), 'nearest', 'extrap');
            end
            % Transform GYR data to match EEG amplitude (for plotting purposes)
            amp_gyro = std(gyrRes,1);
            amp_eeg = mean(std(EEG.data,[],2));
            d = (amp_eeg ./ amp_gyro) / 5;
            gyrRes = gyrRes .* d;
            disp(['GYR amplitudes were multiplied by ' num2str(d(1)) ' (X), ' num2str(d(2)) ' (Y), ' num2str(d(3)) ' (Z) to match EEG data scale.']);
            nChans = size(EEG.data,1);
            EEG.data(nChans+1:nChans+3,:) = gyrRes';
            EEG.chanlocs(nChans+1).labels = 'GYR_X';
            EEG.chanlocs(nChans+2).labels = 'GYR_Y';
            EEG.chanlocs(nChans+3).labels = 'GYR_Z';
            EEG.nbchan = EEG.nbchan+3;
            EEG = eeg_checkset(EEG);
        end
    else
        ind_gyro = find(startsWith(varNames, 'gyro_'));
        ind_gyro = ind_gyro - 1;
        gyroData = table2array(data(:,ind_gyro));

        if strcmp(rec_type, 'muse_monitor')  %same sampling rate as EEG
            gyroData(nans,:) = [];
        else                                    %muse_direct: different sampling rate
            disp('Filling empty values of GYR data to match EEG sampling rate...');
            vals = find(~isnan(gyroData(:,1)));
            for i = 1:size(vals,1)
                Vals = gyroData(vals(i),:);
                if i == 1
                    for j = 1:vals(i)
                        gyroData(j,:) = Vals;
                    end
                else
                    for j = vals(i-1)+1:vals(i)
                        gyroData(j,:) = Vals;
                    end
                end
                if i == size(vals,1)
                    for j = vals(i):size(gyroData,1)
                        gyroData(j,:) = Vals;
                    end
                end
            end
            gyroData(nans,:) = [];    %Remove NaNs to match EEG data length
        end

        %Transform data to match EEG amplitude
        amp_gyro = std(gyroData,1);
        amp_eeg = mean(std(EEG.data,[],2));
        for i = 1:length(amp_gyro)
            d(i) = (amp_eeg / amp_gyro(i)) / 5; %#ok<AGROW>
            gyroData(:,i) = bsxfun(@times, gyroData(:,i), d(i));
        end
        disp(['GYR_X amplitude was multiplied by ' num2str(d(1)) ' to match EEG data scale.']);
        disp(['GYR_Y amplitude was multiplied by ' num2str(d(2)) ' to match EEG data scale.']);
        disp(['GYR_Z amplitude was multiplied by ' num2str(d(3)) ' to match EEG data scale.']);

        %Add to EEG structure
        nChans = size(EEG.data,1);
        EEG.data(nChans+1:nChans+3,:) = gyroData';
        EEG.chanlocs(nChans+1).labels = 'GYR_X';
        EEG.chanlocs(nChans+2).labels = 'GYR_Y';
        EEG.chanlocs(nChans+3).labels = 'GYR_Z';
        EEG.nbchan = EEG.nbchan+3;
        EEG = eeg_checkset(EEG);
    end
end

%% PPG

if params.ppg
    if ~strcmp(rec_type, 'muse_direct')
        warning('PPG import is only implemented for Muse Direct recordings; skipping PPG.');
    else
        disp('Importing PPG data.');
        ind_ppg = find(startsWith(varNames, 'ppg_'));
        ind_ppg = ind_ppg - 1;
        ppgData = table2array(data(:,ind_ppg));

        disp('Resampling PPG data to match EEG sampling rate...');
        vals = find(~isnan(ppgData(:,1)));
        for i = 1:size(vals,1)
            Vals = ppgData(vals(i),:);
            if i == 1
                for j = 1:vals(i)
                    ppgData(j,:) = Vals;
                end
            else
                for j = vals(i-1)+1:vals(i)
                    ppgData(j,:) = Vals;
                end
            end
            if i == size(vals,1)
                for j = vals(i):size(ppgData,1)
                    ppgData(j,:) = Vals;
                end
            end
        end
        ppgData(nans,:) = [];    %Remove NaNs to match EEG data length

        % Correct signal by substracting ambient light from red diode signal
        disp('Removing ambient light signal (PPG1) from blood flow signal (PPG3)');
        ppgData_corr = ppgData(:,3) - ppgData(:,1);

        % Adjust amplitude
        amp_ppg = std(ppgData_corr);
        amp_eeg = mean(std(EEG.data,[],2));
        d = (amp_ppg/amp_eeg)/2;
        ppgData_corr = ppgData_corr ./ d;

        %Add to EEG structure
        nChans = size(EEG.data,1);
        EEG.data(nChans+1,:) = ppgData_corr';   %remove ambient light from PPG signal
        EEG.chanlocs(nChans+1).labels = 'PPG';  %ppg1 = sensor (recording continuously ambient light until diode is ON)
        EEG.nbchan = EEG.nbchan+1;
        EEG = eeg_checkset(EEG);
    end
end

%% AUX

if params.aux
    if ~strcmp(rec_type, 'muse_monitor')
        error('AUX import is only available for MindMonitor (Muse 1) recordings; this recording App does not include AUX data');
    else
        disp('Importing AUX data.');
        % user AUX electrodes are AUX_E/aux_right...; AUX_RIGHT is the Muse reference pin, not an electrode
        ind_aux = find(startsWith(varNames, 'aux') & ~strcmpi(varNames, 'aux_right'));
        if ~isempty(ind_aux)
            ind_aux = ind_aux - 1;
            auxData = table2array(data(:,ind_aux));
            auxData(nans,:) = [];

            %Adjust amplitude
            amp_aux = std(auxData);
            amp_eeg = mean(std(EEG.data,[],2));
            d = amp_aux/amp_eeg;
            auxData = auxData ./ d;

            nChans = size(EEG.data,1);
            if size(auxData,2) == 1
                EEG.data(nChans+1,:) = auxData';
                EEG.chanlocs(nChans+1).labels = 'AUX';
                EEG.nbchan = EEG.nbchan+1;
            else
                if size(auxData,2) > 2
                    warning('More than 2 AUX columns found; importing the first 2.');
                    auxData = auxData(:,1:2);
                end
                EEG.data(nChans+1:nChans+2,:) = auxData';
                EEG.chanlocs(nChans+1).labels = 'AUX1';
                EEG.chanlocs(nChans+2).labels = 'AUX2';
                EEG.nbchan = EEG.nbchan+2;
            end
            EEG = eeg_checkset(EEG);
        else
            warning('No AUX electrodes were used in this recording (no AUX column found); skipping AUX.');
        end
    end
end

%% Resample at the manufacturer default sample rate if needed
if exist('EEG','var') && EEG.srate ~= 256
    if strcmp(rec_type, 'muse_monitor')
        disp('Make sure the sampling rate is set to "Constant" in the settings of your MindMonitor App!');
    end
    warning('Forcing resampling at 256 Hz (i.e., Manufacturer''s default sampling rate)');
    EEG = pop_resample(EEG, 256);
end

disp('MUSE data were imported into EEGLAB.');

%% Detect bad channels using trained classifiers
% Models trained on 5-s windows of raw data filtered 1-50 Hz; this must be
% done on the fly here (on the 4 EEG channels only), or accuracy may drop.

detectBadChan = any(strcmpi(varargin, 'detectBadChan'));

if detectBadChan
    disp('Scanning file to detect bad channels...')
    if EEG.nbchan < 4
        warning('Bad-channel detection requires the 4 Muse EEG channels; skipping.');
    else
        maxTol = .5;    % max portion of bad 5-s windows to tolerate before flagging a channel
        vis = 0;        % set to 1 inside scan_channels call to visualize flagged channels
        TMPEEG = pop_select(EEG, 'channel', 1:4);   % classifiers need only the 4 EEG channels
        [badChan, badChanLabels] = scan_channels(TMPEEG, maxTol, vis);
        if any(badChan)
            EEG = pop_select(EEG,'nochannel',badChanLabels);
        end
    end
end

%% Command history
if nargin < 1
    flagNames = {'acc' 'gyr' 'ppg' 'aux' 'optics'};
    flags = [params.acc params.gyr params.ppg params.aux params.optics];
    optFlags = flagNames(flags == 1);
else
    optFlags = varargin;
end
if isempty(optFlags)
    com = sprintf('EEG = import_muse(''%s'');', file_path);
else
    com = sprintf('EEG = import_muse(''%s'', %s);', file_path, vararg2str(optFlags));
end

end

%% Local subfunctions

% Sample rate as the mode of sample counts per second of recording
function sRate = local_sRateFromTime(Time, nRows)
if isnumeric(Time) && any(abs(Time(1)) > 1e8)
    % numeric POSIX seconds passed directly (Athena packet/OSC formats)
    Time = datetime(Time, 'ConvertFrom', 'posixtime');
end
tNum = datenum(Time);
span = tNum(end) - tNum(1);     % in days; knownRates in Hz
% Preferred estimate: (n-1)/span (robust to timestamps truncated by the app)
% snapped to the nearest known Muse rate when close; else mode of samples/second
knownRates = [64 128 256 512];
if span*86400 > 2 && abs(tNum(end) - tNum(1))*86400 < 1e6 % guard: >= 2 s of data for the span estimate
    if span > 0
        rowRate = (nRows-1) / (span*86400);
        [~, k] = min(abs(rowRate - knownRates));
        if abs(rowRate - knownRates(k)) / knownRates(k) < 0.10
            sRate = knownRates(k);
            disp(['Sample rate detected: ' num2str(sRate) ' Hz (from ' num2str(nRows) ' samples over ' num2str(span*86400, '%.2f') ' s)']);
            return;
        end
    end
    secBins = floor((tNum - tNum(1))*86400) + 1;
    perSec = accumarray(secBins(:), 1);
    perSec = perSec(perSec > 0);
    sRate = mode(perSec);
    disp(['Sample rate detected: ' num2str(sRate) ' Hz (mode over ' num2str(numel(perSec)) ' s)']);
    return;
end
% shorter than 1 s of data: estimate from the total span (timestamps may be
% truncated to ms/whole s by the recording app, which breaks median-dt)
if span <= 0
    error('Could not determine the sample rate (inconsistent timestamps)');
end
sRate = (nRows-1) / (span*86400);
% snap to the nearest known Muse rate when close (timestamps truncation bias)
[~, k] = min(abs(sRate - knownRates));
if abs(sRate - knownRates(k)) / knownRates(k) < 0.20
    sRate = knownRates(k);
else
    sRate = round(sRate);
end
disp(['Sample rate detected: ' num2str(sRate) ' Hz (recording shorter than 1 s)']);
end

% Parse payload rows ('790.45,712.7,718.8,804.5' comma lists) into a numeric matrix
function [M, tSecOut] = local_parseRows(payload, tSec, rows)
nRows = numel(rows);
tSecOut = zeros(nRows, 1);
cellPayload = cell(nRows, 1);
maxVals = 0;
for i = 1:nRows
    p = payload(rows(i));
    tSecOut(i) = tSec(rows(i));
    if isnumeric(p) && isscalar(p)
        cellPayload{i} = {num2str(p)};
    elseif iscellstr(p) || ischar(p) || isstring(p)
        if ischar(p) || isstring(p)
            cellPayload{i} = strsplit(char(string(p)), ',');
        else
            % single cellstr row
            cellPayload{i} = strsplit(char(string(p{1})), ',');
        end
    elseif iscell(p)
        % readtable can return the quoted payload split in multiple cells; concatenate
        joined = strjoin(arrayfun(@(s) char(string(s)), p(:).', 'uni', 0), '');
        cellPayload{i} = strsplit(joined, ',');
    else
        cellPayload{i} = {};
    end
    maxVals = max(maxVals, numel(cellPayload{i}));
end
M = nan(nRows, maxVals);
for i = 1:nRows
    v = str2double(cellPayload{i});
    M(i, 1:numel(v)) = v;
end
end

% Optics channel labels (Mind Monitor mapping for the Muse S Athena)
function labels = local_opticLabels(n)
switch n
    case 4
        labels = {'LI_730' 'RI_730' 'LI_850' 'RI_850'};
    case 8
        labels = {'LO_730' 'RO_730' 'LO_850' 'RO_850' 'LI_730' 'RI_730' 'LI_850' 'RI_850'};
    case 16
        labels = {'LO_730' 'RO_730' 'LO_850' 'RO_850' 'LI_730' 'RI_730' 'LI_850' 'RI_850' ...
                  'LO_Red' 'RO_Red' 'LO_Amb' 'RO_Amb' 'LI_Red' 'RI_Red' 'LI_Amb' 'RI_Amb'};
    otherwise
        labels = arrayfun(@(k) sprintf('Opt%d', k), 1:n, 'uni', 0);
end
end