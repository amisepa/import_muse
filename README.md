# import_muse EEGLAB plugin

This plugin imports Muse .csv files recorded with either the Mind Monitor App or the Muse Direct App. Compatible with Muse 1 (2014 and 2016), Muse 2, and Muse S. Automatically converts data to the EEGLAB format. 

Non-EEG channels (Accelerometer, Gyroscope, Photoplethysmogram, and Auxiliary) can be exported with the EEG data (resampled and slightly transformed to fit), or as separate outputs (raw, untouched). 

This plugin automatically converts each data type to the EEGLAB format, providing access to EEGLAB's advanced tools (e.g. filtering, clean_rawdata, LIMO).

If this plugin does not work for you, see also this other independent implementation for [importing Muse data](https://github.com/sccn/eeglab_musemonitor_plugin).

## Graphic interface

![image](https://user-images.githubusercontent.com/58382227/120024250-bb6d2980-bfa3-11eb-9980-6f6b1b87161f.png)

![](https://github.com/amisepa/import_muse/blob/main/wiki/img30.png)

![](https://github.com/amisepa/import_muse/blob/main/wiki/img35.png)


## Usage

See Wiki (https://github.com/amisepa/import_muse/wiki) for usage and examples.

## Flag bad channels using trained classifers

< EEG = import_muse(filepath,'detectBadChan'); >

The EEG signals in input must be raw (no prior preprocessing), and will be band-pass filtered by this function for best classification performance (but your EEG file output remains raw). The only parameter to change (maxTol) is how much of a channel should be 
tolerated as bad before it is flagged as bad.
For each window, some features are computed (RMS, SNR, low-frequency power), which were selected as most important by a Random Forest model during model training and validation. 

Various ML models were trained and tested (decision trees, logistic regression, LDA, SVM, Naive Bayes, neural networks). They implement PCA-dimension reduction, hyperparameter tuning, and 5-fold cross-validation. Training was done on 80% of a dataset. After model validation, models were tested on the remaining 20% of data (different individuals). The best models achieved 93.5% for frontal channels (logistic regression) and 91.4% for the posterior channels (decision tree).

## Version history
v1.1 - added trained classifiers to flag bad channels
v1.0 - Plugin created and available - June 7, 2021

## Interaxon's Muse specs

Manufacturer website: https://choosemuse.com/muse-2/

![](https://github.com/amisepa/import_muse/blob/main/wiki/img27.png)

![](https://github.com/amisepa/import_muse/blob/main/wiki/img28.png)

![](https://github.com/amisepa/import_muse/blob/main/wiki/img29.png)

## Signal validation (literature)

![](https://github.com/amisepa/import_muse/blob/main/wiki/img31.png)

![](https://github.com/amisepa/import_muse/blob/main/wiki/img32.png)

![](https://github.com/amisepa/import_muse/blob/main/wiki/img33.png)

From Krigolson et al., 2017 (https://doi.org/10.3389/fnins.2017.00109): 
![](https://github.com/amisepa/import_muse/blob/main/wiki/img34.png)

Our validation for frequency domain, peak alpha frequency, and alpha asymmetry: 

Cannard, C., Wahbeh, H., & Delorme, A. (2021, December). Validating the wearable MUSE headset for EEG spectral analysis and Frontal Alpha Asymmetry. In 2021 IEEE International Conference on Bioinformatics and Biomedicine (BIBM) (pp. 3603-3610). IEEE.

https://www.biorxiv.org/content/10.1101/2021.11.02.466989v1.full.pdf
