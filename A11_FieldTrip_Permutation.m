%% Permutation analysis comparing two conditions (baseline to post, or sham to active)

clear variables
clc

addpath /Users/uqcbrad2/Documents/MATLAB/eeglab13_6_5b;
addpath(genpath('/Users/uqcbrad2/Documents/MATLAB/custom_functions'));

% Define paths for working directory
path_general = '/Users/uqcbrad2/Desktop/EEGTest';
path_raw = '/Users/uqcbrad2/Desktop/EEGTest/raw_preproc';
path_ica = '/Users/uqcbrad2/Desktop/EEGTest/interpolated_ica';
path_prod = '/Users/uqcbrad2/Desktop/EEGTest/production';

addpath /Users/uqcbrad2/Documents/MATLAB/fieldtrip-master
ft_defaults

cd(path_ica)% go to where the ica files are
fileList = dir('*denf_i_icasub_nofilt_interval*.set'); % select all the files

%!!!!!!!!!!!!
re_ref='no'; % yes or no
new_ref='avg'; % 'Cz' or 'avg'
%!!!!!!!!!!!!

if strcmp(re_ref,'no')==1
    % Load any file so that we have the electrode names and positions
    eeglab;
    [~, base_filename, ~] = fileparts(fileList(1).name);
    EEG = pop_loadset('filename',[base_filename '.set'],'filepath',path_ica);
    [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, 0 );
    
    cd(path_prod);
    load('tPowerOriginal.mat');
    disp('Original reference used');
    
    %     % Introduce some fake data to check the methods
    %     tPower.A.meanfreq(1:18,1:2,1:5)=100000;
    %     tPower.S.meanfreq(1:1,1:3,1:5)=1;
    
elseif strcmp(re_ref,'yes')==1 && strcmp(new_ref,'avg')==1
    % Load any file so that we have the electrode names and positions
    eeglab;
    [~, base_filename, ~] = fileparts(fileList(1).name);
    EEG = pop_loadset('filename',[base_filename '_ravg.set'],'filepath',path_ica);
    [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, 0 );
    
    cd(path_prod);
    load('tPoweravg.mat');
    disp('Average reference used');
    
elseif strcmp(re_ref,'yes')==1 && strcmp(new_ref,'Cz')==1
    % Load any file so that we have the electrode names and positions
    eeglab;
    [~, base_filename, ~] = fileparts(fileList(1).name);
    EEG = pop_loadset('filename',[base_filename '_rCz.set'],'filepath',path_ica);
    [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, 0 );
    
    cd(path_prod);
    load('tPowerCz.mat');
    disp('Cz reference used');
    
end % end if re_ref

disp('Excluded subjects:'), tPower.excluded_sub

%% Getting the Data structures ready from previous analysis
%!!!!!!!!!!!!!!!!!!!!!!
transit_Data1 = permute (tPower.A.meanfreq,[1 3 2]); % here, you can enter the two conditions you want to compare, for example comparing the average of all post-sham spectra to the baseline sham spectra:tPower.S.meanfreq and tPower.S1.meanfreq; or if comparing mean of post-active and mean of post-sham (in dB): tPower.A.meanfreq_dB and tPower.S.meanfreq_PC
idxToRemove = all(all(isnan(transit_Data1),3),2); % getting rid of participants marked as NaN
transit_Data1(idxToRemove,:,:) = [];

transit_Data2 = permute (tPower.A1.meanfreq,[1 3 2]);
idxToRemove = all(all(isnan(transit_Data2),3),2);
transit_Data2(idxToRemove,:,:) = [];
%!!!!!!!!!!!!!!!!!!!!!!

%% Possible to select subsets of electrodes

chan_to_plot_labels = {EEG.chanlocs.labels};%{'AF3';'AF7';'C2';'C4';'C5';'C6';'CP2';'CP4';'CP5';'CP6';'F1';'F2';'F3';'F4';'F5';'F6';'F7';'F8';'FC2';'FC4';'FC5';'FC6';'FT7';'FT8';'O1';'O2';'Oz';'P1';'P2';'P3';'P4';'P5';'P6';'P7';'P8';'PO3';'PO4';'PO7';'PO8';'POz';'Pz';'TP7';'TP8'};%{EEG.chanlocs.labels};
% whole scalp limited: {'F7';'F3';'Fz';'F4';'F8';'T8';'C4';'T7';'P7';'P3';'Pz';'P4';'P8';'O2';'Oz';'O1'};
% bordering the stimulation electrodes: {'F8';'F6';'F4';'F2';'Fz';'F1';'F3';'F5';'FC5';'C5';'CP5';'CP3';'CP1';'CPz';'FCz'};

chan_to_plot = []; % get the indices for 'plotchans'
for chan_i = 1:length(chan_to_plot_labels)
    idx = find(ismember({EEG.chanlocs.labels},(chan_to_plot_labels(chan_i)))); % find the index for each channel
    if ~isnan(idx)
        chan_to_plot(chan_i) = idx;
    end % end if
end %chan_i

%% Initialise structures

Data1 = []; % Condition1
Data2 = []; % Condition2
X = {EEG.chanlocs.X};
Y = {EEG.chanlocs.Y};
Z = {EEG.chanlocs.Z};
posX=[];posY=[];posZ=[];

% select the data from our electrodes of interest and get their positions
for chan_i=1:length(chan_to_plot)
    Data1(:,chan_i,:) = transit_Data1(:,chan_to_plot(chan_i),:);
    Data2(:,chan_i,:) = transit_Data2(:,chan_to_plot(chan_i),:);
    posX(chan_i) = X{chan_to_plot(chan_i)};
    posY(chan_i) = Y{chan_to_plot(chan_i)};
    posZ(chan_i) = Z{chan_to_plot(chan_i)};
end %chan_i

%% Configuration

% format data
Condition1.label        = chan_to_plot_labels;%{EEG.chanlocs.labels};
Condition1.dimord       = 'subj_chan_freq';
Condition1.freq         = [0.75 2.5 6 10 13.5 20]; % using the center of each frequency band

Condition1.powspctrm    = Data1;
Condition2              = Condition1;
Condition2.powspctrm    = Data2;

% define electrode locations, layout and neighbourhood
sens.label    = chan_to_plot_labels';
sens.elecpos  = [posX; posY; posZ]';
sens.chanpos  = [posX; posY; posZ]'; % we need chanpos in case we wanted to recombine data at different location than sensors
elec = ft_datatype_sens(sens);

cfg         = [];
cfg.elec    = elec;
cfg.rotate  = 90; % somehow the positions of electrodes were rotated
layout      = ft_prepare_layout(cfg);
% ft_layoutplot(cfg); % visualising

cfg                 = [];
cfg.method          = 'distance'; %'distance', 'triangulation' or 'template'; 'triangulation' did not work well for me
cfg.neighbourdist   = 0.18; % maximum distance between neighbouring sensors (only for 'distance'); empirically defined here on our 64 channels
cfg.layout          = layout;
cfg.feedback        ='yes';
myneighbours        = ft_prepare_neighbours(cfg);

%% Perform cluster based permutation analysis

cfg             = [];
cfg.method      = 'montecarlo';
cfg.numrandomization = 1000;
cfg.correctm    = 'cluster'; %'cluster''no', 'max', 'bonferroni', 'holm', 'hochberg', 'fdr'; no or max seem to return the most consistent results (with "fake" data)
cfg.alpha       = 0.025; % for two-tailed
cfg.clusteralpha = 0.05;
cfg.neighbours  = myneighbours;
cfg.design      = [1:size(Data1,1) 1:size(Data1,1); ones(1,size(Data1,1)) 2*ones(1,size(Data1,1))];
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.uvar        = 1;
cfg.ivar        = 2;
% cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details

sigZeroStats_clust = ft_freqstatistics(cfg, Condition1, Condition2);
save('sigZeroStats_clustA_A1.mat','sigZeroStats_clust');

% Plot the result
cfg = [];
cfg.alpha  = 0.05;
cfg.zlim   = [-4 4];
cfg.layout = layout;
ft_clusterplot(cfg, sigZeroStats_clust);
colorbar