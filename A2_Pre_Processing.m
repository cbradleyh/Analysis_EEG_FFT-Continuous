%% Prior to starting this script, create 'StartEndAll.xlsx'

% That file should contain all the start and end times of the segments we want to analyse. 
% These are determined 'by eye' under BrainVision software.
% For now, the file contains at least 9 columns,starting with the filename
% (e.g.01_XXXXX_Group-A_Ses-1_0003) and the start and end times (in
% seconds) of the 4 intervals of interest (before tACS, and after each of the 3 blocks of tACS).

%% Reading files, epoching, downsampling and filtering

clear; close all; clc;

addpath /Users/uqcbrad2/Documents/MATLAB/eeglab13_6_5b;

% Define paths for working directory
path_general = '/Users/uqcbrad2/Desktop/EEGTest';
path_raw = '/Users/uqcbrad2/Desktop/EEGTest/raw_preproc'; % I often create 5 sub-folders to the general folder: ...
% ... bin (for binary files), code (for these scripts and others), interpolated_ica (for the treated files), production (for the resulting graphics, power values, etc...) and raw_prepoc (for the preprocessing steps).

% Fetch file-name and time-intervals to be analysed from an external file
cd(path_general) % go where the xlsx file is
[~,~,interval_data] = xlsread('StartEndAll.xlsx'); % 
cd(path_raw) % go where the raw files are

% Define new sampling rate that will be used if necessary (if there have been different sampling rates in the study, or if you only need 1000Hz for analysis)
new_sample_rate = 1000; % in Hz

%  Set a counting variable for determining subject number
subject_number_count=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, treat one file at a time, by going through the .xlsx spreadsheet one
% line at a time
for file_i = 2:length(interval_data) % we start at 2 because of the header line
    
    % Print on screen which file we are analysing
    fprintf('%s',['Analyzing File ' char(interval_data(file_i,1))]);
    
    % Determine (from the code and session number) whether the stimulation was active or sham
    % !!!!! This will have to change depending on your naming practice !!!!!!
    if ~isempty(strfind(char(interval_data(file_i,1)),'Group-A_Ses-2')) || ~isempty(strfind(char(interval_data(file_i,1)),'Group-B_Ses-1'))% if filename contains Group-A_Ses-2 or Group-B_Ses-1 == active
        tACS_type={'_tACS_Active'};
    elseif ~isempty(strfind(char(interval_data(file_i,1)),'Group-B_Ses-2')) || ~isempty(strfind(char(interval_data(file_i,1)),'Group-A_Ses-1'))
        tACS_type={'_tACS_Sham'};
    end
      
    % Read and define file names and intervals
    filename = char(strcat(interval_data(file_i,1), '.vhdr'));
    new_file_name = char(strcat(interval_data(file_i,1), tACS_type)); % adding the 'active' or 'sham' information to the file name
    
    % !!!!! If you don't have 4 segments, change here !!!!!
    intervals = [interval_data{file_i,2} interval_data{file_i,3}; interval_data{file_i,4} interval_data{file_i,5}; interval_data{file_i,6} interval_data{file_i,7};interval_data{file_i,8} interval_data{file_i,9}]; % select the start and end of the 4 segments (before tACS, and after each of the 3 blocks of tACS)
    
    % Load the dataset to work on
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; % open EEGLAB
    EEG = pop_loadbv(path_raw, filename);
    newchans = convertlocs(EEG.chanlocs, 'sph2all'); % converts Matlab spherical coordinates to all (cartesian 3D, topo, etc...)
    EEG.chanlocs = newchans;
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',new_file_name,'gui','off');
    
    % Downsample EEG data to a reasonable rate
    if EEG.srate > new_sample_rate
        EEG = pop_resample(EEG, new_sample_rate);
    end
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',strcat(new_file_name,'_d'),'gui','off'); % 'd' for 'downsampled'
    
    % Select the intervals you are interested in
        EEG = pop_select( EEG,'time',intervals ); 
        % Apply a notch filter (50Hz) and a low-pass 100Hz filter
        EEG = pop_eegfiltnew(EEG, 49, 51, [], 1, [], 1); % pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder, revfilt, usefft, plotfreqz, minphase)
        EEG = pop_eegfiltnew(EEG, [], 100, [], [], [], 1);
        % Demean
        EEG = pop_rmbase(EEG, [], []); % it helps with visualisation and with exclusion of 'flat' channels later on
        % Save the resulting segments by creating a new set
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',strcat(new_file_name,'_denf'),'savenew',strcat(new_file_name,'_denf.set'),'gui','off'); % 'denf' for 'downsampled, epoched, notched and fitltered'
 
    % just making sure the EEG datasets are not mixed together for the next round
    clear('ALLEEG','EEG','CURRENTSET')
    
end % end file_i

% fclose(fid);
fprintf('%s', 'Done!');
