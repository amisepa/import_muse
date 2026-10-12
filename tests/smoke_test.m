% smoke_test() - Import every recording in test_files/ with import_muse and
% check that the output is a consistent EEGLAB structure. Run by the GitHub
% Actions workflow (.github/workflows/test.yml) and usable locally:
%
%   >> cd import_muse/tests; smoke_test
%
% EEGLAB must be on the MATLAB path, or its folder given in the EEGLAB_PATH
% environment variable.

function smoke_test

pluginPath = fileparts(fileparts(mfilename('fullpath')));
if ~exist('eeglab', 'file')
    eeglabPath = getenv('EEGLAB_PATH');
    if isempty(eeglabPath)
        error('EEGLAB not found: add it to the path or set EEGLAB_PATH');
    end
    addpath(eeglabPath);
end
eeglab('nogui');
addpath(pluginPath, '-begin');   % test this checkout, not an installed copy of the plugin

testPath = fullfile(pluginPath, 'test_files');

% file name, import options
tests = {
    'muse-monitor3.csv'             {}
    'muse-monitor3.csv'             {'acc' 'gyr' 'aux'}
    'muse-direct1.csv'              {}
    'muse-direct1.csv'              {'acc' 'gyr' 'ppg'}
    'muse-direct2.csv'              {}
    'muse-direct2.csv'              {'detectBadChan'}
    'muse-direct2.csv'              {'reref'}
    'athena_eeg_export.csv'         {}
    'athena_musedirect_packet.csv'  {}
    'athena_musedirect_packet.csv'  {'optics'}
    'athena_muselab_osc.csv'        {}
    'athena_muselab_osc.csv'        {'optics'}
    'athena_optics_export.csv'      {'optics'}
    };

for iTest = 1:size(tests,1)
    fileName = tests{iTest,1};
    opts = tests{iTest,2};
    fprintf('\n=== %s %s\n', fileName, strjoin(opts, ' '));

    [EEG, com] = import_muse(fullfile(testPath, fileName), opts{:});
    EEG = eeg_checkset(EEG);

    assert(~isempty(EEG.data), 'no data imported');
    assert(EEG.nbchan == size(EEG.data,1), 'nbchan does not match data');
    assert(EEG.pnts == size(EEG.data,2), 'pnts does not match data');
    assert(EEG.srate > 0, 'invalid sampling rate');
    assert(numel(EEG.chanlocs) == EEG.nbchan, 'chanlocs do not match channels');
    assert(all(isfinite(EEG.data(:))), 'non-finite values in data');
    assert(ischar(com) && ~isempty(com), 'empty history string');

    % options must add channels to the EEG-only import of the same file
    if isempty(opts)
        nChanEEG = EEG.nbchan;
    elseif any(ismember(opts, {'acc' 'gyr' 'aux' 'ppg'}))
        assert(EEG.nbchan > nChanEEG, 'options did not add channels');
    end
end

fprintf('\nAll %d imports passed.\n', size(tests,1));
