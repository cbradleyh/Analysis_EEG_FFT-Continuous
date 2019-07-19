%% Prior to starting this script, you should know which channels have the worst tACS artefacts 

% This should have been done during the visual inspection phase earlier.
% These channels are often the ones on or just around the tACS pads. You
% will enter them below

%% Exclude tACS-artefacted electrodes and interpolate bad channels

clear; close all; clc;

addpath /Users/uqcbrad2/Documents/MATLAB/eeglab13_6_5b;

% Define paths for working directory
path_general = '/Users/uqcbrad2/Desktop/EEGTest';
path_raw = '/Users/uqcbrad2/Desktop/EEGTest/raw_preproc';
path_ica = '/Users/uqcbrad2/Desktop/EEGTest/interpolated_ica';

cd(path_raw) % go to where the raw files are 
fileList = dir('*denf.set'); % select all the files

% Now, treat one file at a time
for file_i = 1:length(fileList)
    
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; % open EEGLAB
    
    % Extract file name information
    [~, base_filename, ~] = fileparts(fileList(file_i).name);
    filename = char(fileList(file_i).name);
    new_filename = [base_filename '_i.set']; % 'i' for 'interpolated'
    
    % Load the dataset to work on
    EEG = pop_loadset(filename,path_raw);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % Exclude tACS-artefacted and systematically noisy electrodes
    % electrodes (TP9 and TP10 and Fpz and Fp1)
    % !!!!!!! This will change depending on your montage and the electrodes
    % you noted are 'contaminated' !!!!!!!!!!
    EEG = pop_select(EEG, 'nochannel', {'FC3' 'FC1' 'C1' 'C3' 'CP1' 'CP3' 'Fp2' 'AF4' 'AF8' 'Fpz' 'Fp1' 'TP9' 'TP10'});
       
    % Save the original EEG locations (minus excluded electrodes) for use in interpolation later
    EEG.allchan = EEG.chanlocs;

    % Identify and interpolate flat electrodes
    flat_chan=[];
    for elec_i = 1:EEG.nbchan
        size_eeg = size(EEG.data,2);
        line = EEG.data(elec_i,size_eeg-1);
        if mean(line)==0
            flat_chan = [flat_chan elec_i];
        end
    end
    EEG = pop_interp(EEG, flat_chan, 'spherical');
    
    % Reject and interpolate bad channels
    % either automatically identified (but check that satisfactory), or
    % from prior visual inspection
    
    % automatic version:
    [EEG, indelec] = pop_rejchan(EEG, 'elec',[1:size(EEG.data,1)] ,'threshold',2,'norm','on','measure','spec','freqrange',[0.2 50]); % these thresholds were determined empirically by checking whether it marked the electrodes of interest in several different recordings.
    EEG.badelec=indelec; % store the information in this field
    EEG = pop_interp(EEG, EEG.allchan, 'spherical'); % interpolate these electrodes

%     % file-based version:
%         % Extract bad channel information from the Excel file and interpolate
%     bad_channels=[];
%     for col_i = 2:size(bad_channels_data,2)
%         if isnan(cell2mat(bad_channels_data(file_i,col_i)));
%             bad_channels_data(file_i,col_i)  = {'NaN'}; % transforming in string otherwise ismember function does not like the input
%         end
%         idx = find(ismember({EEG.chanlocs.labels},(bad_channels_data(file_i,col_i)))); % find the index for each channel
%         if ~isnan(idx)
%             bad_channels(col_i-1) = idx;
%         end % end if
%     end % end col_i
%     EEG = pop_interp(EEG, bad_channels(bad_channels>0), 'spherical');
    
    % Save file
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = pop_saveset( EEG, 'filepath',path_ica,'filename',new_filename);
    
    % just making sure the EEG datasets are not mixed together for the next round
    clear('ALLEEG','EEG','CURRENTSET')
    
end % file_i
fprintf('%s', 'Done interpolation!');

%% Run ICA on all channels
clear; close all; clc;
%!!!!!!!!!!!!!!!
addpath /Users/uqcbrad2/Documents/MATLAB/eeglab13_6_5b;

% Define paths for working directory
path_general = '/Users/uqcbrad2/Desktop/EEGTest';
path_raw = '/Users/uqcbrad2/Desktop/EEGTest/raw_preproc';
path_ica = '/Users/uqcbrad2/Desktop/EEGTest/interpolated_ica';

cd(path_ica) % go to where the raw files are 
fileList = dir('*denf_i.set'); % select all the files

% Treat one file at a time
for file_i = 1:length(fileList)
    
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; % open EEGLAB
    
    % Extract file name information
    [~, base_filename, ~] = fileparts(fileList(file_i).name);
    filename = char(fileList(file_i).name);
    new_filename = [base_filename '_2-35_ICA.set'];
    
    % Load the dataset to work on
    EEG = pop_loadset(filename,path_ica);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % Downsample and filter to reduce processing time and increase quality of ICA results
    EEG = pop_resample(EEG, 256);
    EEG = pop_eegfiltnew(EEG, 2, 35, [], 0, [], 1);
    
    % run the ICA
    EEG = pop_runica(EEG,'icatype','binica','extended',1,'pca',rank(EEG.data)); % do a rank-reducing PCA to make sure we don't use the interpolated channels
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    % Save in a new set
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',new_filename,'savenew',new_filename,'gui','off');
    
    % just making sure the EEG datasets are not mixed together for the next round
    clear('ALLEEG','EEG','CURRENTSET')
    
end % file_i
fprintf('%s', 'Done ICA!');

