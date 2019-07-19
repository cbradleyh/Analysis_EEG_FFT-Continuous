    
clear; close all; clc;

addpath /Users/uqcbrad2/Documents/MATLAB/eeglab13_6_5b;

% Define paths for working directory
path_general = '/Users/uqcbrad2/Desktop/EEGTest';
path_raw = '/Users/uqcbrad2/Desktop/EEGTest/raw_preproc';
path_ica = '/Users/uqcbrad2/Desktop/EEGTest/interpolated_ica';

cd(path_ica)
fileList = dir('*_i_icasub_nofilt.set')%dir('*_i_icasub_nofilt.set'); % select all the files


%% Now, treat one file at a time
for file_i = 1:length(fileList)
    
    %Initialize
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; % open EEGLAB
    A=[];B=[];C=[];D=[];E=[];

    % Extract file name information 
    [~, base_filename, ~] = fileparts(fileList(file_i).name);    
    % Load the dataset
    EEG = pop_loadset(fileList(file_i).name,path_ica);
    [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, 1 );
    
    % segments
    intervals=[];
    for i = 1:length(EEG.event)
        if EEG.event(i).urevent==2
           intervals = [intervals, EEG.event(i).latency];
        end;
    end

    intervals = intervals/1000; % bring back to a unit of seconds
    intervals = [intervals, EEG.xmax]; % add the end point for the last segment
    
    % cycle through segments, saving and naming them as you go
    for i = 1:length(intervals)-1
        
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',1,'study',0); % make sure we are retrieving the full dataset that has just received ICA transformation
        EEG = pop_select( EEG,'time',[intervals(i),intervals(i+1)] ); 
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',strcat(base_filename,'_interval',num2str(i)),'savenew',strcat(base_filename,'_interval',num2str(i),'.set'),'gui','off'); %
        
    end
end
