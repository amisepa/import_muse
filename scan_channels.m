%% Scan Muse channels using 5-s sliding windows and use trained models to 
% classify each window as good/bad automatically and quickly. Frontal and 
% posterior channels are considered separately (one fine-tuned classifier 
% for each because their signal characteristics are pretty different). 
% 
% The EEG signals in input must be raw (no prior preprocessing), and will
% be band-pass filtered by this function for best classification performance
% (but your EEG file output remains raw).
% 
% The only parameter to change (maxTol) how much of a channel should be 
% tolerated as bad before it is flagged as bad.
% 
% For each window, some features are computed (RMS, SNR, low-frequency
% power), which were selected as most important by a Random Forest model 
% during model training and validation. 
% 
% Various ML models were trained and tested (decision trees, logistic 
% regression, LDA, SVM, Naive Bayes, neural networks). They implement 
% PCA-dimension reduction, hyperparameter tuning, and 5-fold
% cross-validation. Training was done on 80% of a dataset. After model 
% validation, models were tested on the remaining 20% of data
% (different individuals). The best models achieved 93.5% for frontal 
% channels (logistic regression) and 91.4% for the posterior channels (decision 
% tree).
% 
%
% Inputs:
%   EEG                 - EEG structure (EEGLAB format)
%   maxTol              - max portion of bad windows to tolerate before 
%                           flagging the channel (default = .5 for 50% of 
%                           windows).
%   vis                 - visualize the flagged channels (1) or not (0)
%
% Output:
%   badChan             - a logical array indicating if channels are bad (1) or
%                           ok (0). A channel is considered bad when the portion of
%                           bad sliding segments is > maxBad.
%   badChanLabels       - cell array with the labels of bad channels
%   badSegs             - logical array indicating the sliding windows
%                           classified as bad by the trained models. 
%
% Example:
%   [badChan, badChanLabels] = scan_channels(EEG, .5, 1);
% 
% Cedric Cannard, March 2023

function [badChan, badChanLabels, badSegs] = scan_channels(EEG,maxTol,vis)

pluginPath = fileparts(which('eegplugin_import_muse.m'));
load(fullfile(pluginPath, 'classifier_front.mat'));  % default = Logistic regression
load(fullfile(pluginPath, 'classifier_post.mat'));  % default = Decision tree

% Defaults
badChanLabels = [];
% if isempty(winSize)
winSize = 5;   % should remain 5 seconds (same as for training)
% end % size of sliding windows in s
winSize = winSize * EEG.srate;          % convert to samples
nSeg = floor(EEG.pnts/winSize);
if isempty(maxTol)
    maxTol = .5; 
end
if isempty(vis)
    vis = false; 
end

% Filter EEG signal (same as for training of classifiers)
EEG = pop_eegfiltnew(EEG,'locutoff',1,'hicutoff',50);    

% Design filter for SNR feature
b = design_fir(100,[2*[0 45 50]/EEG.srate 1],[1 1 0 0]);

% Extract features
badChan = false(1,EEG.nbchan);
disp('Lookg for bad channels using trained models...')
for iChan = 1:EEG.nbchan
    signal = EEG.data(iChan,:);

    flatSeg = nan(nSeg-1,1);
    feat1 = nan(nSeg-1,1);
    feat2 = nan(nSeg-1,1);
    feat3 = nan(nSeg-1,1);
    feat4 = nan(nSeg-1,1);
    feat5 = nan(nSeg-1,1);
    for iWind = 1:nSeg-1

        % Sliding window
        if iWind == 1
            tStart = 1;
        else
            tStart = tEnd + 1;
        end
        tEnd = (tStart + winSize)-1;
        if tEnd > EEG.pnts
            d = tEnd - EEG.pnts;
            tEnd = EEG.pnts;
            warning('tEnd is %g s beyond the last sample. Replacing with last sample.', round(d*256,1))
        end
        tSeg = tStart:tEnd;  % time index in samples for this segment
        
        %%% Frontal features (see trainedModelFront.RequiredVariables)
        
        % RMS
        feat1(iWind,:) = rms(signal(tSeg));

        % LF power
        tmp = get_psd(signal(tSeg),256,'hamming',50,[],256,[0.2 3],'psd');
        feat2(iWind,:) = rms(tmp);
    
        % SNR
        tmp = filtfilt_fast(b,1,signal(tSeg)');
        feat3(iWind,:) = mad(signal(tSeg) - tmp');
        
        % POST features (see trainedModelPost.RequiredVariables)
        if iChan == 1 || iChan == 4 
            
            % Kurtosis
            feat4(iWind,:) = kurtosis(signal(tSeg));

            % HF power
            tmp = get_psd(signal(tSeg),256,'hamming',50,[],256,[70 100],'psd');
            feat5(iWind,:) = rms(tmp);
        end

        % Additional classification: if more than 50% of segment is flat
        flatPortion = sum( abs(diff(signal(tSeg)))<(20*eps) ) / (length(tSeg)-1);
        if flatPortion > .5
            flatSeg(iWind,:) = true;  % 2 for bad
        else
            flatSeg(iWind,:) = false;  % 1 for good
        end
    end
    
    if iChan == 2 || iChan == 3 
        features = table(feat1,feat2,feat3,'variablenames',{ 'RMS' 'LF_power' 'SNR'});
        prediction = trainedModelFront.predictFcn(features); % classify windows
    else
        features = table(feat1,feat2,feat3,feat4,feat5,'variablenames',{ 'RMS' 'LF_power' 'SNR' 'Kurtosis' 'HF_power'});
        prediction = trainedModelPost.predictFcn(features); % classify windows
    end

    % Tag if good or bad
    badSegs(iChan,:) = or(single(prediction) == 2, flatSeg);
    
    % Bad channels
    if sum(badSegs(iChan,:))/(nSeg-1) > maxTol
        badChan(iChan) = true;
    end
end

if sum(badChan) == 0 
    disp('No bad channels detected by trained classifier')
else
    badChanLabels = {EEG.chanlocs(badChan).labels};
    fprintf('Bad channels detected by trained classifier: '); fprintf('%s ', badChanLabels{:}); fprintf('\n')
end

if vis
    % eegplot(EEG.data,'winlength',30,'srate',EEG.srate,'spacing',1000);
    TMPEEG = pop_select(EEG,'nochannel', badChanLabels);
    TMPEEG.etc.clean_channel_mask = ~badChan;
    vis_artifacts(TMPEEG,EEG);
end


%% Subfunctions

% FIR filter design from Christian's clean_artifacts code
function B = design_fir(N,F,A,nfft,W)
if nargin < 4 || isempty(nfft)
    nfft = max(512,2^ceil(log(N)/log(2))); 
end
if nargin < 5
    W = 0.54 - 0.46*cos(2*pi*(0:N)/N); 
end
F = interp1(round(F*nfft),A,(0:nfft),'pchip');
F = F .* exp(-(0.5*N)*sqrt(-1)*pi*(0:nfft)./nfft);
B = real(ifft([F conj(F(end-1:-1:2))]));
B = B(1:N+1).*W(:)';

% filtfilt_fast from Christian's clean_artifacts code
function X = filtfilt_fast(varargin)
% if nargin == 3
[B, A, X] = deal(varargin{:});
% elseif nargin == 4
%     [N, F, M, X] = deal(varargin{:});
%     B = design_fir(N,F,sqrt(M)); A = 1;
% end
if A == 1
    was_single = strcmp(class(X),'single');
    w = length(B); t = size(X,1);    
    X = double([bsxfun(@minus,2*X(1,:),X(1+mod(((w+1):-1:2)-1,t),:)); X; bsxfun(@minus,2*X(t,:),X(1+mod(((t-1):-1:(t-w))-1,t),:))]);
    X = filter_fast(B,A,X); X = X(length(X):-1:1,:);
    X = filter_fast(B,A,X); X = X(length(X):-1:1,:);
    X([1:w t+w+(1:w)],:) = [];
    if was_single
        X = single(X); end    
else    
    X = filtfilt(B,A,X);
end

% filter_fast from Christian's clean_artifacts code
function [X,Zf] = filter_fast(B,A,X,Zi,dim)
if nargin <= 4
    dim = find(size(X)~=1,1); 
end
if nargin <= 3
    Zi = []; 
end
lenx = size(X,dim);
lenb = length(B);
if lenx == 0
    Zf = Zi;
elseif lenb < 256 || lenx<1024 || lenx <= lenb || lenx*lenb < 4000000 || ~isequal(A,1)
    if nargout > 1
        [X,Zf] = filter(B,A,X,Zi,dim);
    else
        X = filter(B,A,X,Zi,dim);
    end
else
    was_single = strcmp(class(X),'single');
    if isempty(Zi)
        if nargout < 2
            X = unflip(oct_fftfilt(B,flip(double(X),dim)),dim);
        else
            X = flip(X,dim);
            [dummy,Zf] = filter(B,1,X(end-length(B)+1:end,:),Zi,1); %#ok<ASGLU>
            X = oct_fftfilt(B,double(X));
            X = unflip(X,dim);
        end
    else
        X = flip(X,dim);
        tmp = filter(B,1,X(1:length(B),:),Zi,1);
        if nargout > 1
            [dummy,Zf] = filter(B,1,X(end-length(B)+1:end,:),Zi,1); %#ok<ASGLU>
        end
        X = oct_fftfilt(B,double(X));
        X(1:length(B),:) = tmp;
        X = unflip(X,dim);
    end
    if was_single
        X = single(X); 
    end
end

function X = flip(X,dim)
if dim ~= 1
    order = 1:ndims(X);
    order = order([dim 1]);
    X = permute(X,order);
end

function X = unflip(X,dim)
if dim ~= 1
    order = 1:ndims(X);
    order = order([dim 1]);
    X = ipermute(X,order);
end


% Compute power spectral density (PSD)
function [pxx, f] = get_psd(eegData,winSize,taperM,overlap,nfft,Fs,fRange,type)
fh = str2func(taperM);
overlap = winSize/(100/overlap); % convert overlap to samples
for iChan = 1:size(eegData,1)
    [pxx(iChan,:), f] = pwelch(eegData(iChan,:),fh(winSize),overlap,nfft,Fs,type);
end
freq = dsearchn(f,fRange(1)):dsearchn(f, fRange(2)); % extract frequencies of interest
f = f(freq(2:end))';
pxx = pxx(:,freq(2:end));     
pxx = 10*log10(pxx); % normalize to deciBels (dB) 
