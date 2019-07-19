%% Apply the ICA results to the original (non-filtered) data-set

clear; close all; clc;

addpath /Users/uqcbrad2/Documents/MATLAB/eeglab13_6_5b;

% Define paths for working directory
path_general = '/Users/uqcbrad2/Desktop/EEGTest';
path_raw = '/Users/uqcbrad2/Desktop/EEGTest/raw_preproc';
path_ica = '/Users/uqcbrad2/Desktop/EEGTest/interpolated_ica';

cd(path_raw) % go to where the raw files are 
fileList = dir('*denf.set'); % select all the files
cd(path_ica)

% Now, treat one file at a time
for file_i = 2:length(fileList)
    
    %Initialize
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; % open EEGLAB
    A=[];B=[];C=[];D=[];E=[];

    % Extract file name information 
    [~, base_filename, ~] = fileparts(fileList(file_i).name);    
    old_filename = [base_filename '_i.set']; % call the set file
    ica_filename = [base_filename '_i_2-35_ICA.set']; % call the set file
    new_filename = [base_filename '_i_icasub_nofilt.set']; % call the set file
    
    % Load the filtered dataset
    EEG = pop_loadset(ica_filename,path_ica);
    [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, 1 );
    % grab the transformation that we have calculated previously
    A=EEG.icawinv;
    B=EEG.icaweights;
    C=EEG.icasphere;
    D=EEG.icachansind;
    E=EEG.reject;
    
    % Load the unfiltered dataset
    EEG = pop_loadset(old_filename,path_ica);
    [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, 2 );
    % Actually apply a high-pass filter at 0.1Hz to get rid of really low
    % variations
    %EEG = pop_eegfiltnew(EEG, 0.1, [], [], [], [], 1);

%     % Make sure that no previously existing ICA information remains
%     EEG.icawinv=[];
%     EEG.icaweights=[];
%     EEG.icasphere=[];
%     EEG.icachansind=[];
%     EEG.reject=[];
%     [ALLEEG, EEG] = eeg_store(ALLEEG, EEG,2);
%     EEG = pop_saveset(EEG,'check', 'on','savemode', 'resave');

    % Replace by new ones
    EEG.icawinv=A;
    EEG.icaweights=B;
    EEG.icasphere=C;
    EEG.icachansind=D;
    EEG.reject=E;
    
    % Subtract components
    EEG = pop_subcomp(EEG,[],0);
    
    % Make sure the boundary events for the segments are well tagged (not
    % tagging the removal of 'noise'
    for i = 1:length(EEG.event)
        EEG.event(i).urevent=2;
    end
        
    % and save again
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, 3);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',new_filename, 'filepath', path_ica);
    
end

fprintf('%s', 'Done ICA correction!');
