# import_muse EEGLAB plugin

This plugin imports Muse .csv files recorded with either the Mind Monitor App or the Muse Direct App, and is compatible with data recorded with the Muse 1 (2014 and 2016), Muse 2, and Muse S.

This plugin automatically converts each Muse data type to the EEGLAB format for easy access to its advanced preprocessing and statistical tools (e.g. filtering, clean_rawdata, LIMO).

Additionally, the plugin allows users to import the following non-EEG channels: Accelerometer (ACC), Gyroscope (GYR), Photoplethysmogram (PPG), or Auxiliary (AUX). 

Users can choose which data channels they wish to import and if they wish to export it with the EEGLAB EEG structure (e.g. data were resampled when sampling rates were different across channels), or as separate outputs (in EEGLAB structure format) for separate analyses. 

# Graphic interface

![image](https://user-images.githubusercontent.com/58382227/120024250-bb6d2980-bfa3-11eb-9980-6f6b1b87161f.png)


# Command line usage

See Wiki (https://github.com/amisepa/import_muse/wiki) for usage and examples.

# Version history
v1.0 - Plugin created and available - June 7, 2021
