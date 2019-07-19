%% This will allow reformating of the large matrix you just created. You can cut it up in different frequency bands, and exclude participants if relevant.
% the current script will need adapting if you don't have 4 segments
% (including a baseline segment)

clear variables
clc

addpath /Users/uqcbrad2/Documents/MATLAB/eeglab13_6_5b;
addpath(genpath('/Users/uqcbrad2/Documents/MATLAB/custom_functions'));

% Define paths for working directory
path_general = '/Users/uqcbrad2/Desktop/EEGTest';
path_raw = '/Users/uqcbrad2/Desktop/EEGTest/raw_preproc';
path_ica = '/Users/uqcbrad2/Desktop/EEGTest/interpolated_ica';
path_prod = '/Users/uqcbrad2/Desktop/EEGTest/production';

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
    load('PowerOriginal.mat');
    
elseif strcmp(re_ref,'yes')==1 && strcmp(new_ref,'avg')==1
    % Load any file so that we have the electrode names and positions
    eeglab;
    [~, base_filename, ~] = fileparts(fileList(1).name);
    EEG = pop_loadset('filename',[base_filename '_ravg.set'],'filepath',path_ica);
    [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, 0 );
    
    cd(path_prod);
    load('Poweravg.mat');
    
elseif strcmp(re_ref,'yes')==1 && strcmp(new_ref,'Cz')==1
    % Load any file so that we have the electrode names and positions
    eeglab;
    [~, base_filename, ~] = fileparts(fileList(1).name);
    EEG = pop_loadset('filename',[base_filename '_rCz.set'],'filepath',path_ica);
    [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, 0 );
    
    cd(path_prod);
    load('PowerCz.mat');
    
end % end if re_ref

%% Define frequency band values and subjects to exclude

% !!!!!!!!!!
Band_values = [0.50,1;1,4;4,8;8,12;12,15;15,25]; % min and max (columns) for each frequency band of interest (rows)
Band_labels = {'tACS frequency', 'delta', 'theta', 'low-alpha', 'high-alpha', 'beta'};
tPower.band.labels = Band_labels;
tPower.band.values = Band_values;

% !!!!!!!!!!
Sub_to_exclude = [2,4,11,12,13,18,21]; % on general data quality from AR phase
tPower.excluded_sub = Sub_to_exclude;

% get a few useful numbers
% !!!!!!!!!!
Num_el = length(Power.channels);
tPower.channels = Power.channels; % prepare a structure to 'trim'
Num_files = length(Power.filenames);
Num_sub = Num_files/8; % There are 8 conditions
Num_sub_trimmed = Num_sub - length(Sub_to_exclude);
Num_band = size(Band_values,1);

%% Reshape the data into matrices corresponding to each condition (8)

% we need to separate the active and sham, per each interval
countingA1=1; % this counting variable allows us to have a smaller number of rows in the new matrices than Num_files
countingA2=1;
countingA3=1;
countingA4=1;
countingS1=1;
countingS2=1;
countingS3=1;
countingS4=1;

for file_i = 1:Num_files
    if isempty(strfind(Power.filenames{file_i},'interval1')) == 0 && isempty(strfind(Power.filenames{file_i},'Active')) == 0 % select each condition using a double negative statement
        Power.A1.filenames{countingA1} = Power.filenames{file_i}; % copy the corresponding data
        Power.A1.spectra(countingA1,:,:) = Power.spectra(file_i,:,:);
        Power.A1.freqs(countingA1,:,:) = Power.freqs(file_i,:,:);
        countingA1=countingA1+1;
        
    elseif isempty(strfind(Power.filenames{file_i},'interval1')) == 0 && isempty(strfind(Power.filenames{file_i},'Sham')) == 0
        Power.S1.filenames{countingS1} = Power.filenames{file_i};
        Power.S1.spectra(countingS1,:,:) = Power.spectra(file_i,:,:);
        Power.S1.freqs(countingS1,:,:) = Power.freqs(file_i,:,:);
        countingS1=countingS1+1;
        
    elseif isempty(strfind(Power.filenames{file_i},'interval2')) == 0 && isempty(strfind(Power.filenames{file_i},'Active')) == 0
        Power.A2.filenames{countingA2} = Power.filenames{file_i};
        Power.A2.spectra(countingA2,:,:) = Power.spectra(file_i,:,:);
        Power.A2.freqs(countingA2,:,:) = Power.freqs(file_i,:,:);
        countingA2=countingA2+1;
        
    elseif isempty(strfind(Power.filenames{file_i},'interval2')) == 0 && isempty(strfind(Power.filenames{file_i},'Sham')) == 0
        Power.S2.filenames{countingS2} = Power.filenames{file_i};
        Power.S2.spectra(countingS2,:,:) = Power.spectra(file_i,:,:);
        Power.S2.freqs(countingS2,:,:) = Power.freqs(file_i,:,:);
        countingS2=countingS2+1;
        
    elseif isempty(strfind(Power.filenames{file_i},'interval3')) == 0 && isempty(strfind(Power.filenames{file_i},'Active')) == 0
        Power.A3.filenames{countingA3} = Power.filenames{file_i};
        Power.A3.spectra(countingA3,:,:) = Power.spectra(file_i,:,:);
        Power.A3.freqs(countingA3,:,:) = Power.freqs(file_i,:,:);
        countingA3=countingA3+1;
        
    elseif isempty(strfind(Power.filenames{file_i},'interval3')) == 0 && isempty(strfind(Power.filenames{file_i},'Sham')) == 0
        Power.S3.filenames{countingS3} = Power.filenames{file_i};
        Power.S3.spectra(countingS3,:,:) = Power.spectra(file_i,:,:);
        Power.S3.freqs(countingS3,:,:) = Power.freqs(file_i,:,:);
        countingS3=countingS3+1;
        
    elseif isempty(strfind(Power.filenames{file_i},'interval4')) == 0 && isempty(strfind(Power.filenames{file_i},'Active')) == 0
        Power.A4.filenames{countingA4} = Power.filenames{file_i};
        Power.A4.spectra(countingA4,:,:) = Power.spectra(file_i,:,:);
        Power.A4.freqs(countingA4,:,:) = Power.freqs(file_i,:,:);
        countingA4=countingA4+1;
        
    elseif isempty(strfind(Power.filenames{file_i},'interval4')) == 0 && isempty(strfind(Power.filenames{file_i},'Sham')) == 0
        Power.S4.filenames{countingS4} = Power.filenames{file_i};
        Power.S4.spectra(countingS4,:,:) = Power.spectra(file_i,:,:);
        Power.S4.freqs(countingS4,:,:) = Power.freqs(file_i,:,:);
        countingS4=countingS4+1;
        
    end % end interval
    
end % end for file_i

%% Exclude participants

% Create the "trimmed" matrices that don't contain subjects to exclude by
% assigning them a value of NaN
tPower.A1 = Power.A1;
tPower.A2 = Power.A2;
tPower.A3 = Power.A3;
tPower.A4 = Power.A4;
tPower.S1 = Power.S1;
tPower.S2 = Power.S2;
tPower.S3 = Power.S3;
tPower.S4 = Power.S4;
for sub_i = 1:length(Sub_to_exclude)
    
    tPower.A1.spectra(Sub_to_exclude(sub_i),:,:) = NaN;
    tPower.A2.spectra(Sub_to_exclude(sub_i),:,:) = NaN;
    tPower.A3.spectra(Sub_to_exclude(sub_i),:,:) = NaN;
    tPower.A4.spectra(Sub_to_exclude(sub_i),:,:) = NaN;
    tPower.S1.spectra(Sub_to_exclude(sub_i),:,:) = NaN;
    tPower.S2.spectra(Sub_to_exclude(sub_i),:,:) = NaN;
    tPower.S3.spectra(Sub_to_exclude(sub_i),:,:) = NaN;
    tPower.S4.spectra(Sub_to_exclude(sub_i),:,:) = NaN;
    
end %sub_i

%% Create average spectra for all three post time-points
% for memory: tPower(sub, band/freq, electrode)

tPower.A = [];
tPower.S = [];

tPower.A.freqs=tPower.A1.freqs;
tPower.S.freqs=tPower.S1.freqs;

% Calculate an average spectrum across all three post time-points
for sub_i = 1:Num_sub
    for chan_i = 1:Num_el
        for freqs_i= 1:length(tPower.A1.spectra)
            tPower.A.spectra(sub_i,freqs_i,chan_i) = (tPower.A2.spectra(sub_i,freqs_i,chan_i)+tPower.A3.spectra(sub_i,freqs_i,chan_i)+tPower.A4.spectra(sub_i,freqs_i,chan_i))/3;
            tPower.S.spectra(sub_i,freqs_i,chan_i) = (tPower.S2.spectra(sub_i,freqs_i,chan_i)+tPower.S3.spectra(sub_i,freqs_i,chan_i)+tPower.S4.spectra(sub_i,freqs_i,chan_i))/3;
        end
    end
end

%% Decibel transformation

% Here calculate a decibel transform of activity: 10*log10(post/pre). 
% We are then looking at an X-fold increase or decrease in values compared to baseline.

for sub_i = 1:Num_sub
    for chan_i = 1:Num_el
        for freqs_i= 1:length(tPower.A1.spectra)
            tPower.A.spectra_dB(sub_i,freqs_i,chan_i) = 10*log10(tPower.A.spectra(sub_i,freqs_i,chan_i)/tPower.A1.spectra(sub_i,freqs_i,chan_i));
            tPower.S.spectra_dB(sub_i,freqs_i,chan_i) = 10*log10(tPower.S.spectra(sub_i,freqs_i,chan_i)/tPower.S1.spectra(sub_i,freqs_i,chan_i));
            
            tPower.A2.spectra_dB(sub_i,freqs_i,chan_i) = 10*log10(tPower.A2.spectra(sub_i,freqs_i,chan_i)/tPower.A1.spectra(sub_i,freqs_i,chan_i));
            tPower.S2.spectra_dB(sub_i,freqs_i,chan_i) = 10*log10(tPower.S2.spectra(sub_i,freqs_i,chan_i)/tPower.S1.spectra(sub_i,freqs_i,chan_i));
            tPower.A3.spectra_dB(sub_i,freqs_i,chan_i) = 10*log10(tPower.A3.spectra(sub_i,freqs_i,chan_i)/tPower.A1.spectra(sub_i,freqs_i,chan_i));
            tPower.S3.spectra_dB(sub_i,freqs_i,chan_i) = 10*log10(tPower.S3.spectra(sub_i,freqs_i,chan_i)/tPower.S1.spectra(sub_i,freqs_i,chan_i));
            tPower.A4.spectra_dB(sub_i,freqs_i,chan_i) = 10*log10(tPower.A4.spectra(sub_i,freqs_i,chan_i)/tPower.A1.spectra(sub_i,freqs_i,chan_i));
            tPower.S4.spectra_dB(sub_i,freqs_i,chan_i) = 10*log10(tPower.S4.spectra(sub_i,freqs_i,chan_i)/tPower.S1.spectra(sub_i,freqs_i,chan_i));
        end
    end
end

% % Calculate summary statistics
% 
% % Calculate the individual mean frequency-band values from the
% % dB-transformed spectra
% for sub_i = 1:Num_sub
%     for chan_i = 1:Num_el
%         for band_i = 1:Num_band
%             temp = tPower.A.spectra_dB(sub_i,:,chan_i);
%             tPower.A.meanfreq_dB(sub_i,band_i,chan_i) = mean(temp(tPower.A.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.A.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
%             
%             temp = tPower.A2.spectra_dB(sub_i,:,chan_i);
%             tPower.A2.meanfreq_dB(sub_i,band_i,chan_i) = mean(temp(tPower.A2.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.A2.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
%             
%             temp = tPower.A3.spectra_dB(sub_i,:,chan_i);
%             tPower.A3.meanfreq_dB(sub_i,band_i,chan_i) = mean(temp(tPower.A3.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.A3.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
%             
%             temp = tPower.A4.spectra_dB(sub_i,:,chan_i);
%             tPower.A4.meanfreq_dB(sub_i,band_i,chan_i) = mean(temp(tPower.A4.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.A4.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
%             
%             temp = tPower.S.spectra_dB(sub_i,:,chan_i);
%             tPower.S.meanfreq_dB(sub_i,band_i,chan_i) = mean(temp(tPower.S.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.S.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
%             
%             temp = tPower.S2.spectra_dB(sub_i,:,chan_i);
%             tPower.S2.meanfreq_dB(sub_i,band_i,chan_i) = mean(temp(tPower.S2.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.S2.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
%             
%             temp = tPower.S3.spectra_dB(sub_i,:,chan_i);
%             tPower.S3.meanfreq_dB(sub_i,band_i,chan_i) = mean(temp(tPower.S3.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.S3.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
%             
%             temp = tPower.S4.spectra_dB(sub_i,:,chan_i);
%             tPower.S4.meanfreq_dB(sub_i,band_i,chan_i) = mean(temp(tPower.S4.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.S4.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
%         end % band_i
%     end % chan_i
% end % sub_i
% 
% % Calculate the group mean band values to feed topoplot
% 
% for band_i = 1:Num_band
%     for chan_i =1:Num_el
%         tPower.A.meanmeanfreq_dB(1,band_i,chan_i) = mean(tPower.A.meanfreq_dB(:,band_i,chan_i),'omitnan');
%         tPower.S.meanmeanfreq_dB(1,band_i,chan_i) = mean(tPower.S.meanfreq_dB(:,band_i,chan_i),'omitnan');
%    end %chan_i
% end % num_band

% Calculate the individual mean frequency-band values from the raw data and
% then the dB transform
for sub_i = 1:Num_sub
    for chan_i = 1:Num_el
        for band_i = 1:Num_band
            temp = tPower.A.spectra(sub_i,:,chan_i);
            tPower.A.meanfreq(sub_i,band_i,chan_i) = mean(temp(tPower.A.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.A.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
            
            temp = tPower.A1.spectra(sub_i,:,chan_i);
            tPower.A1.meanfreq(sub_i,band_i,chan_i) = mean(temp(tPower.A1.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.A1.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
            tPower.A.meanfreq_dB(sub_i,band_i,chan_i)= 10*log10(tPower.A.meanfreq(sub_i,band_i,chan_i)/tPower.A1.meanfreq(sub_i,band_i,chan_i));
            
            temp = tPower.A2.spectra(sub_i,:,chan_i);
            tPower.A2.meanfreq(sub_i,band_i,chan_i) = mean(temp(tPower.A2.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.A2.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
            tPower.A2.meanfreq_dB(sub_i,band_i,chan_i)= 10*log10(tPower.A2.meanfreq(sub_i,band_i,chan_i)/tPower.A1.meanfreq(sub_i,band_i,chan_i));
            
            temp = tPower.A3.spectra(sub_i,:,chan_i);
            tPower.A3.meanfreq(sub_i,band_i,chan_i) = mean(temp(tPower.A3.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.A3.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
            tPower.A3.meanfreq_dB(sub_i,band_i,chan_i)= 10*log10(tPower.A3.meanfreq(sub_i,band_i,chan_i)/tPower.A1.meanfreq(sub_i,band_i,chan_i));

            temp = tPower.A4.spectra(sub_i,:,chan_i);
            tPower.A4.meanfreq(sub_i,band_i,chan_i) = mean(temp(tPower.A4.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.A4.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
            tPower.A4.meanfreq_dB(sub_i,band_i,chan_i)= 10*log10(tPower.A4.meanfreq(sub_i,band_i,chan_i)/tPower.A1.meanfreq(sub_i,band_i,chan_i));
                        
            temp = tPower.S.spectra(sub_i,:,chan_i);
            tPower.S.meanfreq(sub_i,band_i,chan_i) = mean(temp(tPower.S.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.S.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
            
            temp = tPower.S1.spectra(sub_i,:,chan_i);
            tPower.S1.meanfreq(sub_i,band_i,chan_i) = mean(temp(tPower.S1.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.S1.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
            tPower.S.meanfreq_dB(sub_i,band_i,chan_i)= 10*log10(tPower.S.meanfreq(sub_i,band_i,chan_i)/tPower.S1.meanfreq(sub_i,band_i,chan_i));
            
            temp = tPower.S2.spectra(sub_i,:,chan_i);
            tPower.S2.meanfreq(sub_i,band_i,chan_i) = mean(temp(tPower.S2.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.S2.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
            tPower.S2.meanfreq_dB(sub_i,band_i,chan_i)= 10*log10(tPower.S2.meanfreq(sub_i,band_i,chan_i)/tPower.S1.meanfreq(sub_i,band_i,chan_i));
            
            temp = tPower.S3.spectra(sub_i,:,chan_i);
            tPower.S3.meanfreq(sub_i,band_i,chan_i) = mean(temp(tPower.S3.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.S3.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
            tPower.S3.meanfreq_dB(sub_i,band_i,chan_i)= 10*log10(tPower.S3.meanfreq(sub_i,band_i,chan_i)/tPower.S1.meanfreq(sub_i,band_i,chan_i));
            
            temp = tPower.S4.spectra(sub_i,:,chan_i);
            tPower.S4.meanfreq(sub_i,band_i,chan_i) = mean(temp(tPower.S4.freqs(sub_i,:,chan_i) >= Band_values(band_i,1) & tPower.S4.freqs(sub_i,:,chan_i) <= Band_values(band_i,2)));
            tPower.S4.meanfreq_dB(sub_i,band_i,chan_i)= 10*log10(tPower.S4.meanfreq(sub_i,band_i,chan_i)/tPower.S1.meanfreq(sub_i,band_i,chan_i));
        
        end % band_i
    end % chan_i
end % sub_i

% Calculate the group mean band values to feed topoplot

for band_i = 1:Num_band
    for chan_i =1:Num_el
        tPower.A1.meanmeanfreq(1,band_i,chan_i) = mean(tPower.A1.meanfreq(:,band_i,chan_i),'omitnan');
        tPower.S1.meanmeanfreq(1,band_i,chan_i) = mean(tPower.S1.meanfreq(:,band_i,chan_i),'omitnan');
        tPower.A.meanmeanfreq(1,band_i,chan_i) = mean(tPower.A.meanfreq(:,band_i,chan_i),'omitnan');
        tPower.S.meanmeanfreq(1,band_i,chan_i) = mean(tPower.S.meanfreq(:,band_i,chan_i),'omitnan');
        tPower.A.meanmeanfreq_dB(1,band_i,chan_i) = mean(tPower.A.meanfreq_dB(:,band_i,chan_i),'omitnan');
        tPower.S.meanmeanfreq_dB(1,band_i,chan_i) = mean(tPower.S.meanfreq_dB(:,band_i,chan_i),'omitnan');
   end %chan_i
end % num_band

%% Saving

if strcmp(re_ref,'no')==1
    save('tPowerOriginal.mat','tPower');
    
elseif strcmp(re_ref,'yes')==1 && strcmp(new_ref,'avg')==1
    save('tPoweravg.mat','tPower');
    
elseif strcmp(re_ref,'yes')==1 && strcmp(new_ref,'Cz')==1
    save('tPowerCz.mat','tPower');
    
end % end if re_ref