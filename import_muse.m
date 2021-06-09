%% import_muse() - Import Muse data recorded with the Muse Monitor or Muse Direct Apps
% Compatible with Muse 1 (models 2014 & 2016), Muse 2, and Muse S
%
% Usage:
%   EEG = import_muse;           % pop-up GUI mode
%   EEG = import_muse(file_path);   %Input filepath as command line and use import/export defaults
%   EEG = import_muse(file_path, import_EEG, import_ACC, import_GYR, import_PPG, import_AUX, export_option);  % command line
%
% Examples:
%   EEG = import_muse(file_path, 1, 1, 1, 0, 1, 2);                          %Import everything and export with EEG (Muse Monitor)
%   EEG = import_muse(file_path, 1, 1, 1, 1, 0, 2);                          %Import everything and export with EEG (Muse Direct)
%   [EEG, ACC, GYR, PPG, AUX] = import_muse(file_path, 1, 1, 1, 0, 1, 1);    %Import everything as separate outputs (Muse Monitor)
%   [EEG, ACC, GYR, PPG, AUX] = import_muse(file_path, 1, 1, 1, 1, 0, 2);    %Import everything as separate outputs (Muse Direct)
%
% Optional input:
%   file_path       - full path and name of the Muse file (e.g. 'path\file_name.csv')
%   import_EEG      - Import EEG data (default ON = 1)
%   import_ACC      - Import accelerometer data (default OFF = 0)
%   import_GYR      - Import gyrometer data (default OFF = 0)
%   import_PPG      - Import photoplethysmogram (PPG) data (default OFF = 0; ONLY works with Muse 2 recorded with Muse Direct).
%   import_AUX      - Import auxiliary channel data (default OFF = 0; for custom-made AUX electrodes)
%   export_option   - Choose to export non-EEG data with the EEG structure (default = 1) or as separate outputs (2)
%
% Outputs:
%   EEG     - EEG data in EEGLAB structure format
%   ACC     - Accelerometer data
%   GYR     - Gyroscope data
%   PPG     - Photoplethysmogram data
%               PPG1 = continuous recording of ambient light until diode is ON
%               PPG2 = IR (SPO2; i.e. oxygen concentration)
%               PPG3 = red diode (blood flow)
%   AUX     - Auxiliary electrode data
%   com     - history string for command line
%
% Author: Cedric Cannard, CerCo, CNRS
%
% Copyright (C) 2021 Cedric Cannard
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version. This program is distributed in
% the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG, ACC, GYR, PPG, AUX, com] = import_muse(file_path, varargin)

com = '';

%File path and inputs
if nargin < 1
    [fileName, filePath] = uigetfile2({'*.csv'}, 'Select Muse .csv file - import_muse()');
    file_path = fullfile(filePath, fileName);
end
if nargin < 2
    export_choice = {'With the EEG data (default)', 'As seperate outputs'};
    uilist = {
        {'style' 'checkbox' 'string' 'Import EEG data (TP9, AF7, AF8, TP10)' 'tag' 'eeg' 'value' 1 'enable' 'on' } ...
        {} ...
        {'Style' 'text' 'string' 'Non-EEG channels: '} ...
        {'style' 'checkbox' 'string' 'Import Accelerometer (ACC) data (x,y,z)' 'tag' 'acc' 'value' 0 'enable' 'on' } ...
        {'style' 'checkbox' 'string' 'Import Gyroscope (GYR) data (x,y,z)' 'tag' 'gyr' 'value' 0 'enable' 'on' } ...
        {'style' 'checkbox' 'string' 'Import Photoplethysmogram (PPG) data (Muse 2 recorded with Muse Direct ONLY)' 'tag' 'ppg' 'value' 0 'enable' 'on' } ...
        {'style' 'checkbox' 'string' 'Import Auxiliary (AUX) data (Muse Monitor only)' 'tag' 'aux' 'value' 0 'enable' 'on' } ...
        {} ...
        {'Style' 'text' 'string' 'Exportation choice for non-EEG channels: '} ...
        {'Style' 'popupmenu' 'String' export_choice 'Tag' 'export'} ...
        };
    uigeom = { 1 1 1 1 1 1 1 1 1 1};
    opt = inputgui(uigeom, uilist, 'pophelp(''import_muse'')', 'Import Muse data -- import_muse()');
    if isempty(opt), return; end
    
    % Defaults
elseif nargin < 7 && nargin > 1
    opt{1} = 1;
    opt{2} = 0;
    opt{3} = 0;
    opt{4} = 0;
    opt{5} = 0;
    opt{6} = 1;
else
    opt = varargin;
end

%% IMPORT DATA AND VARIABLE NAMES

disp('Importing data...');
data = readtable(file_path);
variable = importdata(file_path);
varNames = variable.textdata;
if size(varNames,2) == 1
    varNames = lower(split(varNames(1,1), ','))';
elseif size(varNames,1) > 1 && size(varNames,2) > 1
    varNames = cellstr(lower(split(varNames(1,1), ',')))';
end

%% TIMESTAMPS

ind_times = contains(varNames, 'time');
Time = table2array(data(:,ind_times));
if isdatetime(Time)
    rec_type = 'Muse_Monitor';
    disp('Recording App detected: Muse Monitor');
    if size(data,2) == 2                            %muse_monitor2
        Time = datetime(variable.textdata(2:end), 'Format','HH:mm:ss:SSS');
    else                                            %muse_monitor1
        Time = datetime(Time, 'Format','HH:mm:ss:SSS');
    end
elseif ~isdatetime(Time) && iscell(Time)      %muse_monitor3
    rec_type = 'Muse_Monitor';
    Time = datetime(Time, 'InputFormat', 'mm:ss.S', 'Format','HH:mm:ss:SSS');
elseif ~isdatetime(Time) && isnumeric(Time)   %muse_direct1 and 2
    rec_type = 'Muse_Direct';
    disp('Recording App detected: Muse Direct');
    Time = datetime(Time, 'ConvertFrom','posixtime', 'Format','HH:mm:ss:SSS');
else
    errordlg('Timestamps format not recognized! Please submit an issue on Github and attach your raw file.','Data Error');
end

%% EEG DATA

varNames = varNames(2:end); %remove time column for indices to match data matrices
EEG = [];

if opt{1}
    if strcmp(rec_type, 'Muse_Monitor')     %muse_monitor1, 2, 3
        ind_eeg = contains(varNames, 'raw');
    elseif strcmp(rec_type, 'Muse_Direct')  %muse_direct1
        ind_eeg = contains(varNames, {'eeg'});
        if sum(ind_eeg) > 4
            ind_eeg(find(ind_eeg,1,'first')+4:length(ind_eeg)) = false;
        end
    end
    if size(data,2) == 2    %muse_monitor2
        eegData = variable.data(:,ind_eeg);
    else                    %muse_monitor1, 3, muse_direct2
        data(:,1) = [];
        eegData = table2array(data(:,ind_eeg));
    end
    
    
    %Remove empty rows
    nans = isnan(eegData(:,1));
    eegData(nans,:) = [];
    eegTime = Time;
    eegTime(nans,:) = [];
    
    %Calculate sample rate and test stability
    fprintf('Calculating EEG sampling rate and testing its stability...\n');
    lags = second(eegTime);
    diff_secs = find(logical(diff(round(lags))));
    sRates = diff(diff_secs);
    eeg_sRate = mode(sRates);
    if std(sRates) > 1
        warning(['EEG sample rate unstable and deviates up to ' num2str(round(std(sRates))) ...
            ' samples/s across file. This is likely due to Bluetooth streaming and should not be an issue.'])
    end
    disp(['Sample rate detected: ' num2str(eeg_sRate) ' Hz']);
    if eeg_sRate < 50 && strcmp(rec_type, 'Muse_Monitor')
        warning(['Your sample rate was set to 1 sample/s in your Muse Monitor settings! ' ...
            'It is highly recommended to set it to "Constant" (i.e. 256 Hz).']);
    elseif eeg_sRate > 256
        disp('Setting EEG sample rate to hardware specifications: 256 Hz')
        eeg_sRate = 256;
    end
    
    %Different method to test sample rate stability
    % fprintf('Interpolating sampling rate and testing its stability...\n');
    % nSamples = round(1./seconds(diff(eegTime)));      %time between each sample
    % nSamples(isinf(nSamples)) = 0.001;                %convert inf values to 0.001
    % sRate = round(mode(nSamples));                      % main sRate over file (in Hz)
    % diff_sRates = unique(nSamples);                     %diffrent sample rates across file (in Hz)
    % sRate_weight = round(sum(nSamples == sRate)/length(eegTime)*100,3);
    % if 100 - sRate_weight > 1
    %     warning([num2str(round(100-sRate_weight)) '% of samples have ' num2str(length(diff_sRates)) ' different sample rates! This is a serious problem from recording'])
    %     warning('Using default sample rate of 256 Hz, but your data might be corrupted!');
    %     sRate = 256;    %ADD SPECIFIC SRATE DEPENDING ON MUSE 1 < 2016 or MUSE 1 2016 or Muse 2
    % elseif sRate < 128
    %     warning('Sample rate interpolated from timestamps is too low (<128 Hz)! Using default sample rate of 256 Hz, but your data might be corrupted!');
    %     sRate = 256;    %ADD SPECIFIC SRATE DEPENDING ON MUSE 1 < 2016 or MUSE 1 2016 or Muse 2
    % end
    
    %EEGLAB STRUCTURE
    EEG = eeg_emptyset;
    EEG.chanlocs = struct('labels', {'TP9' 'AF7'   'AF8'   'TP10'});
    % eegData = bsxfun(@minus, eegData, mean(eegData,1));   %substract mean
    EEG.data = eegData';
    EEG.srate = eeg_sRate;
    EEG.pnts   = size(EEG.data,2);
    EEG.nbchan = size(EEG.data,1);
    EEG.xmin = 0;
    EEG.trials = 1;
    EEG.setname = ['EEG data (' rec_type ')'];
    EEG = eeg_checkset(EEG);
    
end

%% ACC

ACC = [];
if opt{2}
    ind_acc = contains(varNames, 'acc');
    if size(data,2) == 2                    %muse_monitor2
        accData = variable.data(:,ind_acc);
    else                                    %muse_monitor1, 3, muse_direct2
        accData = table2array(data(:,ind_acc));
    end
    
    if opt{6} == 2           %Export as a separate output
        ACC = eeg_emptyset;
        ind = ~isnan(accData(:,1));
        ACC.data = accData(ind,:)';
        
        %Interpolate sample rate
        lags = second(Time(ind));
        diff_secs = find(logical(diff(round(lags))));
        sRates = diff(diff_secs);
        sRate = mode(sRates);
        disp(['ACC sample rate detected: ' num2str(sRate)]);
        if sRate > 256
            ACC.srate = 256;
            disp('Setting ACC sample rate to hardware specifications: 256 Hz')
        else
            ACC.srate = sRate;
        end
        ACC.chanlocs = struct('labels', {'ACC_X' 'ACC_Y'   'ACC_Z'});
        ACC.pnts   = size(ACC.data,2);
        ACC.nbchan = size(ACC.data,1);
        ACC.xmin = 0;
        ACC.trials = 1;
        ACC.setname = 'Accelerometer data';
        ACC = eeg_checkset(ACC);
        
    elseif opt{6} == 1       %Export with EEG structure
        if strcmp(rec_type, 'Muse_Monitor')     %same sampling rate
            accData(nans,:) = [];
            
        elseif strcmp(rec_type, 'Muse_Direct') %different sampling rate
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
            accData(nans,:) = [];    %Remove NaNs from eegData to match data length
        end
        
        %Transform ACC data to match EEG amplitude
        amp_acc = mean(accData,1);
        amp_eeg = mean(mean(EEG.data,2));
        for i = 1:length(amp_acc)
            d(i) = (amp_eeg / amp_acc(i)) / 3 ;
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
        % ACC = pop_eegfiltnew(ACC, 'hicutoff', 15);  %for Muse monitor data to smooth signal
    end
end

%% GYRO

GYR = [];
if opt{3}
    ind_gyro = contains(varNames, 'gyro');
    if size(data,2) == 2                    %muse_monitor2
        gyroData = variable.data(:,ind_gyro);
    else                                    %muse_monitor1, 3, muse_direct2
        gyroData = table2array(data(:,ind_gyro));
    end
    
    if opt{6} == 2          %Export as a separate output
        GYR = eeg_emptyset;
        ind = ~isnan(gyroData(:,1));
        GYR.data = gyroData(ind,:)';
        
        %Interpolate sample rate
        lags = second(Time(ind));
        diff_secs = find(logical(diff(round(lags))));
        sRates = diff(diff_secs);
        sRate = mode(sRates);
        disp(['GYR sample rate detected: ' num2str(sRate)]);
        if sRate > 256
            GYR.srate = 256;
            disp('Setting GYR sample rate to hardware specifications: 256 Hz')
        else
            GYR.srate = sRate;
        end
        GYR.chanlocs = struct('labels', {'GYR_X' 'GYR_Y'   'GYR_Z'});
        GYR.pnts   = size(GYR.data,2);
        GYR.nbchan = size(GYR.data,1);
        GYR.xmin = 0;
        GYR.trials = 1;
        GYR.setname = 'Gyroscope data';
        GYR = eeg_checkset(GYR);
        
    elseif opt{6} == 1           %export with EEG structure
        if strcmp(rec_type, 'Muse_Monitor') %same sampling rate
            gyroData(nans,:) = [];
        elseif strcmp(rec_type, 'Muse_Direct') %different sampling rate
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
            gyroData(nans,:) = [];    %Remove NaNs from eegData to match data length
        end
        
        %Transform data to match EEG amplitude
        amp_gyro = mean(gyroData,1);
        amp_eeg = mean(mean(EEG.data,2));
        for i = 1:length(amp_gyro)
            d(i) = (amp_eeg / amp_gyro(i)) / 5;
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

PPG = [];
if opt{4} && strcmp(rec_type, 'Muse_Direct')
    ind_ppg = contains(varNames, 'ppg');
    if size(data,2) == 2                    %muse_monitor2
        ppgData = variable.data(:,ind_ppg);
    else                                    %muse_monitor1, 3, muse_direct2
        ppgData = table2array(data(:,ind_ppg));
    end
    
    if opt{6} == 2       %Export as a separate output
        PPG = eeg_emptyset;
        ind = ~isnan(ppgData(:,1));
        PPG.data = ppgData(ind,:)';
        
        %Interpolate sample rate
        lags = second(Time(ind));
        diff_secs = find(logical(diff(round(lags))));
        sRates = diff(diff_secs);
        sRate = mode(sRates);
        disp(['PPG sample rate detected: ' num2str(sRate)]);
        if sRate > 256
            PPG.srate = 256;
            disp('Setting PPG sample rate to hardware specifications: 256 Hz')
        else
            PPG.srate = sRate;
        end
        PPG.chanlocs = struct('labels', {'PPG1' 'PPG2'   'PPG3'});
        PPG.pnts   = size(PPG.data,2);
        PPG.nbchan = size(PPG.data,1);
        PPG.xmin = 0;
        PPG.trials = 1;
        PPG.setname = 'Photoplethysmogram data';
        PPG = eeg_checkset(PPG);
        
    elseif opt{6} == 1       %Export with EEG structure
        disp('Filling empty values of PPG data to match EEG sampling rate...');
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
        ppgData(nans,:) = [];    %Remove NaNs from eegData to match data length
        
        %Correct signal by substracting ambient light from red diode signal
        disp('Removing ambient light signal (PPG1) from blood flow signal (PP3)');
        ppgData_corr = ppgData(:,3) - ppgData(:,1);
        
        %Adjust amplitude
        amp_ppg = std(ppgData_corr);
        amp_eeg = mean(std(EEG.data,[],2));
        d = (amp_ppg/amp_eeg)/2;
        ppgData_corr = ppgData_corr ./ d;

%         amp_ppg = mean(ppgData,1);
%         amp_eeg = mean(mean(EEG.data,2));
%         for i = 1:length(amp_ppg)
%             %             if i == 1
%             d(i) = amp_ppg(i) / amp_eeg;
%             ppgData(:,i) = ppgData(:,i) ./ d(i);
%             %             else
%             %                 d(i) = (amp_ppg(i) / amp_eeg) / 10;
%             %                 ppgData(:,i) = ppgData(:,i) ./ d(i);
%             %             end
%         end
%         disp(['PPG1 amplitude was divided by ' num2str(d(1)) ' to match EEG data scale.']);
%         disp(['PPG2 amplitude was divided by ' num2str(d(2)) ' to match EEG data scale.']);
%         disp(['PPG3 amplitude was divided by ' num2str(d(3)) ' to match EEG data scale.']);
        
        %Add to EEG structure
        nChans = size(EEG.data,1);
        EEG.data(nChans+1,:) = ppgData_corr';   %remove ambient light from PPG signal
        EEG.chanlocs(nChans+1).labels = 'PPG';     %ppg1 = sensor (recording continuously ambient light until diode is ON)
%         EEG.data(nChans+1:nChans+3,:) = ppgData';
%         EEG.chanlocs(nChans+1).labels = 'PPG1';     %ppg1 = sensor (recording continuously ambient light until diode is ON)
%         EEG.chanlocs(nChans+2).labels = 'PPG2';     %ppg2 = IR (SPO2; i.e. oxygen concentration)
%         EEG.chanlocs(nChans+3).labels = 'PPG3';     %ppg3 = red diode (blood flow)
%         EEG.nbchan = EEG.nbchan+3;
        EEG.nbchan = EEG.nbchan+1;
        EEG = eeg_checkset(EEG);
    end
%     warning('For analysis: PPG1 (i.e. ambient light signal) needs to be substracted from PPG3 (i.e. the red diode signal measuring blood flow) to obtain the correct PPG signal');
    
elseif opt{4} && strcmp(rec_type, 'Muse_Monitor')
    error('Muse Monitor does not record PPG data (let us know if this changes in an update of Muse Monitor)');
end

%% AUX

AUX = [];
if opt{5} && strcmp(rec_type, 'Muse_Monitor')
    ind_aux = contains(varNames, 'aux');
    if size(data,2) == 2                    %muse_monitor2
        auxData = variable.data(:,ind_aux);
    else                                    %muse_monitor1, 3, muse_direct2
        auxData = table2array(data(:,ind_aux));
    end
    auxData(nans,:) = [];
    
    if opt{6} == 2   %Export as a separate output
        AUX = eeg_emptyset;
        AUX.data = auxData';
        AUX.srate = EEG.srate;
        AUX.chanlocs = struct('labels', {'AUX'});
        AUX.pnts   = size(AUX.data,2);
        AUX.nbchan = size(AUX.data,1);
        AUX.xmin = 0;
        AUX.trials = 1;
        AUX.setname = 'Auxiliary data';
        AUX = eeg_checkset(AUX);
        
    elseif opt{6} == 1   %Export with EEG structure
        
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
            EEG.data(nChans+1:nChans+2,:) = auxData';
            EEG.chanlocs(nChans+1).labels = 'AUX1';
            EEG.chanlocs(nChans+1).labels = 'AUX2';
            EEG.nbchan = EEG.nbchan+2;
        end
        EEG = eeg_checkset(EEG);
    end
    
elseif opt{5} && strcmp(rec_type, 'Muse_Direct')
    error('Muse Direct does not record AUX data');
end

%% History

if opt{6} == 1
    com = sprintf('EEG = import_muse(''%s'', %s, %s, %s, %s, %s, %s);', file_path, ...
        vararg2str(opt{1}), vararg2str(opt{2}), vararg2str(opt{3}), vararg2str(opt{4}), vararg2str(opt{5}), vararg2str(opt{6}));
else
    com = sprintf('[EEG, ACC, GYR, PPG, AUX] = import_muse(''%s'', %s, %s, %s, %s, %s, %s);', file_path, ...
        vararg2str(opt{1}), vararg2str(opt{2}), vararg2str(opt{3}), vararg2str(opt{4}), vararg2str(opt{5}), vararg2str(opt{6}));
end

end
