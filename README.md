# import_muse EEGLAB plugin

This plugin imports Muse .csv files recorded with either the Mind Monitor App or the Muse Direct App, and is compatible with data recorded with the Muse 1 (2014 and 2016), Muse 2, and Muse S.

This plugin automatically converts each Muse data type to the EEGLAB format for easy access to its advanced preprocessing and statistical tools (e.g. filtering, clean_rawdata, LIMO).

Additionally, the plugin allows users to import the following non-EEG channels: Accelerometer (ACC), Gyroscope (GYR), Photoplethysmogram (PPG), or Auxiliary (AUX). 

Users can choose which data channels they wish to import and if they wish to export it with the EEGLAB EEG structure (e.g. data were resampled when sampling rates were different across channels), or as separate outputs (in EEGLAB structure format) for separate analyses. 

# Graphic interface

![image](https://user-images.githubusercontent.com/58382227/120024250-bb6d2980-bfa3-11eb-9980-6f6b1b87161f.png)


# Command line usage

eeglab

EEG = import_muse;                                                       %Input everything with GUI

or 

EEG = import_muse(file_path);                                            %Input filepath and select options with GUI)

EEG = import_muse(file_path, 1);                                         %Import EEG only

EEG = import_muse(file_path, 1, 1, 1, 0, 1, 1);                          %Import everything and export with EEG (Muse Monitor)

EEG = import_muse(file_path, 1, 1, 1, 1, 0, 1);                          %Import everything and export with EEG (Muse Direct)

[EEG, ACC, GYR, PPG, AUX] = import_muse(file_path, 1, 1, 1, 0, 1, 2);    %Import everything as separate outputs (Muse Monitor)

[EEG, ACC, GYR, PPG, AUX] = import_muse(file_path, 1, 1, 1, 1, 0, 2);    %Import everything as separate outputs (Muse Direct)


See Wiki (https://github.com/amisepa/import_muse/wiki) for usage and examples.

# Version history
v1.0 - Plugin created and available - June 7, 2021
