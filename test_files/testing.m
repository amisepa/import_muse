% testing import_muse

clear;clc;close all;
eeglab
chanLocs = fileparts(which('dipfitdefs.m'));

file_name = 'muse-direct2.csv';
file_path = fullfile('C:\Users\IONSLAB\Documents\MATLAB\eeglab\plugins\import_muse1.0\test_files', file_name);

% EEG = pop_muse;                                                       %Input everything with GUI
% EEG = pop_muse(file_path);                                            %Input filepath and select options with GUI)
% EEG = pop_muse(file_path, 1, 1, 1, 0, 1, 1);                          %Import everything and export with EEG (Muse Monitor)
% EEG = pop_muse(file_path, 1, 1, 1, 1, 0, 1);                          %Import everything and export with EEG (Muse Direct)
% [EEG, ACC, GYR, PPG, AUX] = pop_muse(file_path, 1, 1, 1, 0, 1, 2);    %Import everything as separate outputs (Muse Monitor)
[EEG, ACC, GYR, PPG, AUX] = import_muse(file_path, 1, 1, 1, 1, 0, 2);    %Import everything as separate outputs (Muse Direct)

% EEG = pop_chanedit(EEG, 'lookup', fullfile(chanLocs, 'standard_BEM', 'elec', 'standard_1020.elc'));
EEG = pop_eegfiltnew(EEG, 'locutoff',1);
% EEG = pop_eegfiltnew(EEG, 'hicutoff',50,'plotfreqz',0);
% figure; pop_spectopo(EEG, 1, [], 'EEG', 'freq', [6 10 22], 'freqrange', [1 50],'electrodes','on');

pop_eegplot(EEG, 1, 1, 1);
pop_eegplot(ACC, 1, 1, 1);
pop_eegplot(GYR, 1, 1, 1);
% pop_eegplot(AUX, 1, 1, 1);
pop_eegplot(PPG, 1, 1, 1);

