%% Do a simple FFT analysis
clear; close all; clc;

addpath /Users/uqcbrad2/Documents/MATLAB/eeglab13_6_5b;

% Define paths for working directory
path_general = '/Users/uqcbrad2/Desktop/EEGTest';
path_raw = '/Users/uqcbrad2/Desktop/EEGTest/raw_preproc';
path_ica = '/Users/uqcbrad2/Desktop/EEGTest/interpolated_ica';
path_prod = '/Users/uqcbrad2/Desktop/EEGTest/production';

cd(path_ica) % go to where the raw files are 
fileList = dir('*denf_i_icasub_nofilt_interval*.set'); % select all the files
cd(path_ica)

%!!!!!!!!!!!!
re_ref='no'; % yes or no
new_ref='avg'; % 'Cz' or 'avg'
%!!!!!!!!!!!!

% Now, treat one file at a time
for file_i = 1:length(fileList)

    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % open EEGLAB
    
    % Extract file name information
    [~, base_filename, ~] = fileparts(fileList(file_i).name);
    filename = char(fileList(file_i).name);
    
    if strcmp(re_ref,'no')==1
        filename = [base_filename '.set']; % call the set file
        EEG = pop_loadset(filename,path_ica); % Load the dataset to work on
        Power.filenames{file_i}= filename;
        Power.channels = {EEG.chanlocs.labels};
        
    elseif strcmp(re_ref,'yes')==1
        filename = [base_filename '.set']; % call the set file
        new_filename = [base_filename '_r' new_ref '.set'];
        EEG = pop_loadset(filename,path_ica); % Load the dataset to work on
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        EEG = eeg_checkset( EEG );
        
        % Re-reference to Cz or average reference, or none
        idx = find(ismember({EEG.chanlocs.labels}, new_ref)); % find indice for the new reference
        EEG = pop_reref(EEG, idx); % idx if we want Cz, [] if average reference
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, 1); % Save file
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset( EEG,'filepath',path_ica,'filename',new_filename);
        Power.filenames{file_i}= new_filename;
        Power.channels = {EEG.chanlocs.labels};
        
    else
        disp 'Do you want to re-reference or not?'
        
    end % end if re_ref
    
    % define some useful variables
    % !!!!!!!!! You may want to change this !!!!!!!!!!!
    srate = EEG.srate; % sampling rate of 1 kHz
    window_length = 5*srate; % Size of the window on which fft will be run, in data points. Initially 10s: 10*srate
    n_overlap = 50*window_length/100; % Overlap between successive windows, in data points, ex: 50%
    
    % Remove the baseline
    EEG.data = rmbase(EEG.data); % in case the ICA transform, filtering and re-referencing changed things - it helps reduce the 'zero-Hz' component
    
    % For each of the selected channels...
    for chan_i = 1:EEG.nbchan
        
        % Calculate the Power Spectrum Density and store it
        [spectra, freqs] = spectopo(EEG.data(chan_i,:),0,srate,'plot','off','winsize',window_length,'overlap',n_overlap);
        convert_spectra = 10.^(spectra/10); % bring spectra back to power/Hz
        Power.spectra(file_i,:,chan_i) = convert_spectra;
        Power.freqs(file_i,:,chan_i) = freqs;
        
    end % chan_i
    
    % just making sure the EEG datasets are not mixed together for the next round
    clear('ALLEEG','EEG','CURRENTSET')
        
end % file_i

% and save information
cd(path_prod) 
if strcmp(re_ref,'no')==1
    save('PowerOriginal.mat','Power');
elseif strcmp(re_ref,'yes')==1
    save(['Power' new_ref '.mat'],'Power')
end

fprintf('%s', 'Done Re-ref and FFT!');