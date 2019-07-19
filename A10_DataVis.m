%% This script allows plotting of data
% Some plotting options are very slow (especially the bar graphs) when
% using all conditions and electrodes. Comment and un-comment to select the
% plotting you want to do.

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

cd(path_prod)

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
    
elseif strcmp(re_ref,'yes')==1 && strcmp(new_ref,'avg')==1
    % Load any file so that we have the electrode names and positions
    eeglab;
    [~, base_filename, ~] = fileparts(fileList(1).name);
    EEG = pop_loadset('filename',[base_filename '_ravg.set'],'filepath',path_ica);
    [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, 0 );
    
    cd(path_prod);
    load('tPoweravg.mat');
    
elseif strcmp(re_ref,'yes')==1 && strcmp(new_ref,'Cz')==1
    % Load any file so that we have the electrode names and positions
    eeglab;
    [~, base_filename, ~] = fileparts(fileList(1).name);
    EEG = pop_loadset('filename',[base_filename '_rCz.set'],'filepath',path_ica);
    [ALLEEG, EEG] = eeg_store( ALLEEG, EEG, 0 );
    
    cd(path_prod);
    load('tPowerCz.mat');
    
end % end if re_ref


%% Possible to select subsets of electrodes

% !!!!!!!!!!!!!
chan_to_plot_labels = {'FCz';'Cz';'CPz';'CP5';'C5';'FC5';'F3';'F1';'P1';'P3'};%{'AF3';'AF7';'C2';'C4';'C5';'C6';'CP2';'CP4';'CP5';'CP6';'F1';'F2';'F3';'F4';'F5';'F6';'F7';'F8';'FC2';'FC4';'FC5';'FC6';'FT7';'FT8';'O1';'O2';'Oz';'P1';'P2';'P3';'P4';'P5';'P6';'P7';'P8';'PO3';'PO4';'PO7';'PO8';'POz';'Pz';'TP7';'TP8'};%{EEG.chanlocs.labels};%{'Fz';'Cz';'Pz'};{'FCz';'Cz';'CPz';'CP5';'C5';'FC5';'F3';'F1';'P1';'P3'};%
% whole scalp limited: {'F7';'F3';'Fz';'F4';'F8';'T8';'C4';'T7';'P7';'P3';'Pz';'P4';'P8';'O2';'Oz';'O1'};
% bordering the stimulation electrodes: {'F8';'F6';'F4';'F2';'Fz';'F1';'F3';'F5';'FC5';'C5';'CP5';'CP3';'CP1';'CPz';'FCz'};

chan_to_plot = []; % get the indices for 'plotchans'
for chan_i = 1:length(chan_to_plot_labels)
    idx = find(ismember({EEG.chanlocs.labels},(chan_to_plot_labels(chan_i)))); % find the index for each channel
        if ~isnan(idx)
            chan_to_plot(chan_i) = idx;
        end % end if
end %chan_i

%% Extract frequency band values and subjects to exclude

Band_values = tPower.band.values; % min and max (columns) for each frequency band of interest (rows)
Band_labels = tPower.band.labels;
Num_el = length(chan_to_plot);
Num_band = length(tPower.band.labels);
% !!!!! Enter how many participants you are supposed to have
Num_sub = 32;
Num_sub_trimmed = Num_sub - length(tPower.excluded_sub);
disp('Excluded subjects:'), tPower.excluded_sub

%% Plot frequency spectra 

% %Compare baseline spectra
% figure
% for chan_i = 1:Num_el
%     X=[];
%     Y=[];
%     EB=[];
%     subplot(ceil(sqrt(length(chan_to_plot))),ceil(sqrt(length(chan_to_plot))),chan_i)
%     X=smooth(mean(tPower.A1.freqs(:,:, chan_to_plot(chan_i)),'omitnan'));
%     Y(1,:)=smooth(mean(10*log10(tPower.A1.spectra(:,:,chan_to_plot(chan_i))),'omitnan'));
%     Y(2,:)=smooth(mean(10*log10(tPower.S1.spectra(:,:,chan_to_plot(chan_i))),'omitnan'));
% 
%     EB(1,:)= smooth(std(10*log10(tPower.A1.spectra(:,:,chan_to_plot(chan_i))),'omitnan')/sqrt(Num_sub_trimmed));
%     EB(2,:)= smooth(std(10*log10(tPower.S1.spectra(:,:,chan_to_plot(chan_i))),'omitnan')/sqrt(Num_sub_trimmed));
%     mseb(X,Y,EB);
%     
%     set(gca,'xlim',[0.5 25],'ylim',[-5 15])
%     title(chan_to_plot_labels{chan_i});
%     
%     yL = get(gca,'YLim');
%     for i=1:length(Band_values)
%         line([Band_values(i,1),Band_values(i,1)],yL,'Color','k');
%     end
%     
% end % end chan_i
% title(['Baseline tACS Active vs Sham ' chan_to_plot_labels{chan_i}]);
% legend({'Active';'Sham'})
% 
% 
% % And now Active to baseline
% figure
% for chan_i = 1:Num_el
%     X=[];
%     Y=[];
%     EB=[];
%     subplot(ceil(sqrt(length(chan_to_plot))),ceil(sqrt(length(chan_to_plot))),chan_i)
%     X=smooth(mean(tPower.A1.freqs(:,:, chan_to_plot(chan_i)),'omitnan'));
%     Y(1,:)=smooth(mean(10*log10(tPower.A.spectra(:,:,chan_to_plot(chan_i))),'omitnan'));
%     Y(2,:)=smooth(mean(10*log10(tPower.A1.spectra(:,:,chan_to_plot(chan_i))),'omitnan'));
% 
%     EB(1,:)= smooth(std(10*log10(tPower.A.spectra(:,:,chan_to_plot(chan_i))),'omitnan')/sqrt(Num_sub_trimmed));
%     EB(2,:)= smooth(std(10*log10(tPower.A1.spectra(:,:,chan_to_plot(chan_i))),'omitnan')/sqrt(Num_sub_trimmed));
%     mseb(X,Y,EB);
%     
%     set(gca,'xlim',[0.5 25],'ylim',[-5 15])
%     title(chan_to_plot_labels{chan_i});
%     
%     yL = get(gca,'YLim');
%     for i=1:length(Band_values)
%         line([Band_values(i,1),Band_values(i,1)],yL,'Color','k');
%     end
%     
% end % end chan_i
% title(['Post tACS Active vs Baseline ' chan_to_plot_labels{chan_i}]);
% legend({'Active';'Baseline'})
%     
% % Zoomed in
% figure
% for chan_i = 1:Num_el
%     X=[];
%     Y=[];
%     EB=[];
%     subplot(ceil(sqrt(length(chan_to_plot))),ceil(sqrt(length(chan_to_plot))),chan_i)
%     X=smooth(mean(tPower.A1.freqs(:,:, chan_to_plot(chan_i)),'omitnan'));
%     Y(1,:)=smooth(mean(10*log10(tPower.A.spectra(:,:,chan_to_plot(chan_i))),'omitnan'));
%     Y(2,:)=smooth(mean(10*log10(tPower.A1.spectra(:,:,chan_to_plot(chan_i))),'omitnan'));
% 
%     EB(1,:)= smooth(std(10*log10(tPower.A.spectra(:,:,chan_to_plot(chan_i))),'omitnan')/sqrt(Num_sub_trimmed));
%     EB(2,:)= smooth(std(10*log10(tPower.A1.spectra(:,:,chan_to_plot(chan_i))),'omitnan')/sqrt(Num_sub_trimmed));
%     mseb(X,Y,EB);
%     
%     set(gca,'xlim',[0.5 5],'ylim',[5 20])
%     title(chan_to_plot_labels{chan_i});
%     
%     yL = get(gca,'YLim');
%     line([1,1],yL,'Color','r');
%     
% end % end chan_i
% title(['Post tACS Active vs Baseline ' chan_to_plot_labels{chan_i}]);
% legend({'Active';'Baseline'})
% 
% % And now Sham to baseline
% figure
% for chan_i = 1:Num_el
%     X=[];
%     Y=[];
%     EB=[];
%     subplot(ceil(sqrt(length(chan_to_plot))),ceil(sqrt(length(chan_to_plot))),chan_i)
%     X=smooth(mean(tPower.A1.freqs(:,:, chan_to_plot(chan_i)),'omitnan'));
%     Y(1,:)=smooth(mean(10*log10(tPower.S.spectra(:,:,chan_to_plot(chan_i))),'omitnan'));
%     Y(2,:)=smooth(mean(10*log10(tPower.S1.spectra(:,:,chan_to_plot(chan_i))),'omitnan'));
% 
%     EB(1,:)= smooth(std(10*log10(tPower.S.spectra(:,:,chan_to_plot(chan_i))),'omitnan')/sqrt(Num_sub_trimmed));
%     EB(2,:)= smooth(std(10*log10(tPower.S1.spectra(:,:,chan_to_plot(chan_i))),'omitnan')/sqrt(Num_sub_trimmed));
%     mseb(X,Y,EB);
%     
%     set(gca,'xlim',[0.5 25],'ylim',[-5 15])
%     title(chan_to_plot_labels{chan_i});
%     
%     yL = get(gca,'YLim');
%     for i=1:length(Band_values)
%         line([Band_values(i,1),Band_values(i,1)],yL,'Color','k');
%     end
%     
% end % end chan_i
% title(['Post tACS Sham vs Baseline ' chan_to_plot_labels{chan_i}]);
% legend({'Sham';'Baseline'})
% 
% % Zoomed in
% figure
% for chan_i = 1:Num_el
%     X=[];
%     Y=[];
%     EB=[];
%     subplot(ceil(sqrt(length(chan_to_plot))),ceil(sqrt(length(chan_to_plot))),chan_i)
%     X=smooth(mean(tPower.S1.freqs(:,:, chan_to_plot(chan_i)),'omitnan'));
%     Y(1,:)=smooth(mean(10*log10(tPower.S.spectra(:,:,chan_to_plot(chan_i))),'omitnan'));
%     Y(2,:)=smooth(mean(10*log10(tPower.S1.spectra(:,:,chan_to_plot(chan_i))),'omitnan'));
% 
%     EB(1,:)= smooth(std(10*log10(tPower.S.spectra(:,:,chan_to_plot(chan_i))),'omitnan')/sqrt(Num_sub_trimmed));
%     EB(2,:)= smooth(std(10*log10(tPower.S1.spectra(:,:,chan_to_plot(chan_i))),'omitnan')/sqrt(Num_sub_trimmed));
%     mseb(X,Y,EB);
%     
%     set(gca,'xlim',[0.5 5],'ylim',[5 20])
%     title(chan_to_plot_labels{chan_i});
%     
%     yL = get(gca,'YLim');
%     line([1,1],yL,'Color','r');
%     
% end % end chan_i
% title(['Post tACS Sham vs Baseline ' chan_to_plot_labels{chan_i}]);
% legend({'Sham';'Baseline'})
% % 
% % Comparing pre-post dB change for both active and sham
%  figure
% for chan_i = 1:Num_el
%     X=[];
%     Y=[];
%     EB=[];
%     subplot(ceil(sqrt(length(chan_to_plot))),ceil(sqrt(length(chan_to_plot))),chan_i)
%     X=smooth(mean(tPower.A1.freqs(:,:, chan_to_plot(chan_i)),'omitnan'));
%     Y(1,:)=smooth(mean(tPower.A.spectra_dB(:,:,chan_to_plot(chan_i)),'omitnan'));
%     Y(2,:)=smooth(mean(tPower.S.spectra_dB(:,:,chan_to_plot(chan_i)),'omitnan'));
%     
%     EB(1,:)= smooth(std(tPower.A.spectra_dB(:,:,chan_to_plot(chan_i)),'omitnan')/sqrt(Num_sub_trimmed));
%     EB(2,:)= smooth(std(tPower.S.spectra_dB(:,:,chan_to_plot(chan_i)),'omitnan')/sqrt(Num_sub_trimmed));
%     mseb(X,Y,EB);
%     
%     set(gca,'xlim',[0.5 5],'ylim',[-1 4]) % yLim=20 with plot
%     title(chan_to_plot_labels{chan_i});
%     
%     yL = get(gca,'YLim');
%     for i=1:length(Band_values)
%         line([Band_values(i,1),Band_values(i,1)],yL,'Color','k');
%     end
%     
% end % end chan_i
% legend({'Active';'Sham'})
% title(['tACS Active vs Sham '  chan_to_plot_labels{chan_i}]);
% % 
%% Mega electrode

% Create and plot mega-electrode

tPower.A1.freqs(:,:,length(EEG.chanlocs)+1)=tPower.A1.freqs(:,:,length(EEG.chanlocs));
tPower.S1.freqs(:,:,length(EEG.chanlocs)+1)=tPower.S1.freqs(:,:,length(EEG.chanlocs));
for sub_i = 1:Num_sub
        for freqs_i= 1:length(tPower.A1.spectra)
            tPower.A.spectra(sub_i,freqs_i,length(EEG.chanlocs)+1) = mean(tPower.A.spectra(sub_i,freqs_i,chan_to_plot),3,'omitnan');
            tPower.A1.spectra(sub_i,freqs_i,length(EEG.chanlocs)+1) = mean(tPower.A1.spectra(sub_i,freqs_i,chan_to_plot),3,'omitnan');
            tPower.S.spectra(sub_i,freqs_i,length(EEG.chanlocs)+1) = mean(tPower.S.spectra(sub_i,freqs_i,chan_to_plot),3,'omitnan');
            tPower.S1.spectra(sub_i,freqs_i,length(EEG.chanlocs)+1) = mean(tPower.S1.spectra(sub_i,freqs_i,chan_to_plot),3,'omitnan');
        end
end

figure

X=[];
Y=[];
EB=[];
X=smooth(mean(tPower.A1.freqs(:,:, length(EEG.chanlocs)+1),'omitnan'));
Y(1,:)=smooth(mean(10*log10(tPower.A.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan'));
Y(2,:)=smooth(mean(10*log10(tPower.A1.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan'));

EB(1,:)= smooth(std(10*log10(tPower.A.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan')/sqrt(Num_sub_trimmed));
EB(2,:)= smooth(std(10*log10(tPower.A1.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan')/sqrt(Num_sub_trimmed));
lineProps.col = {'r','k'};
mseb(X,Y,EB,lineProps);

set(gca,'xlim',[0.5 25],'ylim',[-3 15])
title('Mega electrode');

yL = get(gca,'YLim');
for i=1:length(Band_values)
    line([Band_values(i,1),Band_values(i,1)],yL,'Color','k');
end

title('Post tACS Active vs Baseline');
legend({'Active';'Baseline'})

% and the zoom
figure
X=[];
Y=[];
EB=[];
X=smooth(mean(tPower.A1.freqs(:,:, length(EEG.chanlocs)+1),'omitnan'));
Y(1,:)=smooth(mean(10*log10(tPower.A.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan'));
Y(2,:)=smooth(mean(10*log10(tPower.A1.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan'));

EB(1,:)= smooth(std(10*log10(tPower.A.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan')/sqrt(Num_sub_trimmed));
EB(2,:)= smooth(std(10*log10(tPower.A1.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan')/sqrt(Num_sub_trimmed));
lineProps.col = {'r','k'};
mseb(X,Y,EB,lineProps);

set(gca,'xlim',[0.5 1],'ylim',[10 20])
title('Mega electrode');

yL = get(gca,'YLim');
line([1,1],yL,'Color','r');

title('Post tACS Active vs Baseline ME');
legend({'Active';'Baseline'})

% And now Sham to baseline
figure

X=[];
Y=[];
EB=[];
X=smooth(mean(tPower.A1.freqs(:,:, length(EEG.chanlocs)+1),'omitnan'));
Y(1,:)=smooth(mean(10*log10(tPower.S.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan'));
Y(2,:)=smooth(mean(10*log10(tPower.S1.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan'));

EB(1,:)= smooth(std(10*log10(tPower.S.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan')/sqrt(Num_sub_trimmed));
EB(2,:)= smooth(std(10*log10(tPower.S1.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan')/sqrt(Num_sub_trimmed));
lineProps.col = {'b','k'};
mseb(X,Y,EB,lineProps);

set(gca,'xlim',[0.5 25],'ylim',[-3 15])
title('Mega electrode');

yL = get(gca,'YLim');
for i=1:length(Band_values)
    line([Band_values(i,1),Band_values(i,1)],yL,'Color','k');
end

title('Post tACS Sham vs Baseline');
legend({'Sham';'Baseline'})

% and the zoom
figure
X=[];
Y=[];
EB=[];
X=smooth(mean(tPower.S1.freqs(:,:, length(EEG.chanlocs)+1),'omitnan'));
Y(1,:)=smooth(mean(10*log10(tPower.S.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan'));
Y(2,:)=smooth(mean(10*log10(tPower.S1.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan'));

EB(1,:)= smooth(std(10*log10(tPower.S.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan')/sqrt(Num_sub_trimmed));
EB(2,:)= smooth(std(10*log10(tPower.S1.spectra(:,:,length(EEG.chanlocs)+1)),'omitnan')/sqrt(Num_sub_trimmed));
lineProps.col = {'b','k'};
mseb(X,Y,EB,lineProps);

set(gca,'xlim',[0.5 1],'ylim',[10 20])
title('Mega electrode');

yL = get(gca,'YLim');
line([1,1],yL,'Color','r');

title('Post tACS Sham vs Baseline ME');
legend({'Sham';'Baseline'})

% Bar graphs

% for the mega-electrode
for sub_i = 1:Num_sub
    for band_i= 1:Num_band
        tPower.A1.meanfreq(sub_i,band_i,length(EEG.chanlocs)+1) = mean(tPower.A1.meanfreq(sub_i,band_i,chan_to_plot),3,'omitnan');
        tPower.A2.meanfreq(sub_i,band_i,length(EEG.chanlocs)+1) = mean(tPower.A2.meanfreq(sub_i,band_i,chan_to_plot),3,'omitnan');
        tPower.A3.meanfreq(sub_i,band_i,length(EEG.chanlocs)+1) = mean(tPower.A3.meanfreq(sub_i,band_i,chan_to_plot),3,'omitnan');
        tPower.A4.meanfreq(sub_i,band_i,length(EEG.chanlocs)+1) = mean(tPower.A4.meanfreq(sub_i,band_i,chan_to_plot),3,'omitnan');
        tPower.S1.meanfreq(sub_i,band_i,length(EEG.chanlocs)+1) = mean(tPower.S1.meanfreq(sub_i,band_i,chan_to_plot),3,'omitnan');
        tPower.S2.meanfreq(sub_i,band_i,length(EEG.chanlocs)+1) = mean(tPower.S2.meanfreq(sub_i,band_i,chan_to_plot),3,'omitnan');
        tPower.S3.meanfreq(sub_i,band_i,length(EEG.chanlocs)+1) = mean(tPower.S3.meanfreq(sub_i,band_i,chan_to_plot),3,'omitnan');
        tPower.S4.meanfreq(sub_i,band_i,length(EEG.chanlocs)+1) = mean(tPower.S4.meanfreq(sub_i,band_i,chan_to_plot),3,'omitnan');

    end
end
 
for band_i = 1:Num_band
    figure
        group_mean = [mean(tPower.A1.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan'),...
            mean(tPower.A2.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan'),...
            mean(tPower.A3.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan'),...
            mean(tPower.A4.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan');...
            mean(tPower.S1.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan'),...
            mean(tPower.S2.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan'),...
            mean(tPower.S3.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan'),...
            mean(tPower.S4.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan')];
        group_std = [std(tPower.A1.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan'),...
            std(tPower.A2.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan'),...
            std(tPower.A3.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan'),...
            std(tPower.A4.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan');...
            std(tPower.S1.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan'),...
            std(tPower.S2.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan'),...
            std(tPower.S3.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan'),...
            std(tPower.S4.meanfreq(:,band_i,length(EEG.chanlocs)+1),'omitnan')];
        group_sem = group_std/sqrt(Num_sub_trimmed);

        barwitherr(group_sem,group_mean)
        colormap(gray)
        set(gca,'XTickLabel',{'Active tACS','Sham tACS'})
        ylabel(['Mean Power ' Band_labels{band_i} ' (microV2)'])
        title(['Mean Power ' Band_labels{band_i} ' in mega electrode ']);

    legend({'PretACS';'After 5min';'After 10min';'After 15min'})

end % band_i


%% Diagnose baseline
% % % 
% for band_i=1:Num_band
%     figure
%     for chan_i = 1:length(chan_to_plot)
%                 subplot(ceil(sqrt(length(chan_to_plot))),ceil(sqrt(length(chan_to_plot))),chan_i)
%         
%                 plot(tPower.A1.meanfreq(:,band_i,chan_to_plot(chan_i)),tPower.S1.meanfreq(:,band_i,chan_to_plot(chan_i)),'.')  % tPower(sub, band, electrode)
%                 refline(1,0) % slope 1, intercept 0
%                 labelpoints(tPower.A1.meanfreq(:,band_i,chan_to_plot(chan_i)),tPower.S1.meanfreq(:,band_i,chan_to_plot(chan_i)),[1:Num_sub],'SE', 0.3,1,'outliers_lin', {1,0,'SD',2.5})
%                 %axis equal,       
%                 title([chan_to_plot_labels{chan_i}]);
%             
%     end % end chan_i
%            title([Band_labels(band_i)]);    
% end % band_i

%% Scatter plot per frequency band
% % 
% for band_i=1:Num_band
%     figure
%     for chan_i = 1:length(chan_to_plot)
%                 subplot(ceil(sqrt(length(chan_to_plot))),ceil(sqrt(length(chan_to_plot))),chan_i)
%         
%                 scatter([1:Num_sub], tPower.A1.meanfreq(:,band_i,chan_to_plot(chan_i))) % tPower(sub, band, electrode)
%                 hold on
%                 scatter([1:Num_sub], tPower.S1.meanfreq(:,band_i,chan_to_plot(chan_i)))
%         
%                 title(['Before tACS Active vs Sham ' Band_labels(band_i) chan_to_plot_labels{chan_i}]);
%             
%     end % end chan_i
%                 legend({'Active';'Sham'})
% 
% end % band_i
% 
% %See if dB changes in active relate to sham
% for band_i=1:Num_band
%     figure
%     for chan_i = 1:length(chan_to_plot)
%                 subplot(ceil(sqrt(length(chan_to_plot))),ceil(sqrt(length(chan_to_plot))),chan_i)
%         
%                 plot(tPower.A.meanfreq_dB(:,band_i,chan_to_plot(chan_i)),tPower.S.meanfreq_dB(:,band_i,chan_to_plot(chan_i)),'.')  % tPower(sub, band, electrode)
%                 refline(1,0) % slope 1, intercept 0
%                 labelpoints(tPower.A.meanfreq_dB(:,band_i,chan_to_plot(chan_i)),tPower.S.meanfreq_dB(:,band_i,chan_to_plot(chan_i)),[1:Num_sub],'SE', 0.3,1,'outliers_lin', {1,0,'SD',2.5})
%                 %axis equal,       
%                 title([chan_to_plot_labels{chan_i}]);
%             
%     end % end chan_i
%            title([Band_labels(band_i)]);    
% end % band_i

%notBoxplot();
%% Plot the bar graphs of all conditions
% 
% for band_i = 1:Num_band
%     figure
%     for chan_i = 1:Num_el
%         group_mean = [mean(tPower.A1.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
%             mean(tPower.A2.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
%             mean(tPower.A3.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
%             mean(tPower.A4.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan');...
%             mean(tPower.S1.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
%             mean(tPower.S2.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
%             mean(tPower.S3.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
%             mean(tPower.S4.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan')];
%         group_std = [std(tPower.A1.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
%             std(tPower.A2.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
%             std(tPower.A3.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
%             std(tPower.A4.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan');...
%             std(tPower.S1.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
%             std(tPower.S2.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
%             std(tPower.S3.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
%             std(tPower.S4.meanfreq(:,band_i,chan_to_plot(chan_i)),'omitnan')];
%         group_sem = group_std/sqrt(Num_sub_trimmed);
% 
%         subplot(ceil(sqrt(Num_el)),ceil(sqrt(Num_el)),chan_i)
%         barwitherr(group_sem,group_mean)
%         colormap(gray)
%         set(gca,'XTickLabel',{'Active tACS','Sham tACS'})
%         ylabel(['Mean Power ' Band_labels{band_i} ' (microV2)'])
%         title(['Mean Power ' Band_labels{band_i} ' in electrode ' chan_to_plot_labels{chan_i}]);
%         
%     end % chan_i
%     legend({'PretACS';'After 5min';'After 10min';'After 15min'})
% 
% end % band_i
% % 
% % % Plot the bar graphs of all conditions PercentChange from baseline
% % for band_i = 1:Num_band
% %     figure
% %     for chan_i = 1:Num_el
% %         group_mean = [mean(tPower.A1.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
% %             mean(tPower.A2.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
% %             mean(tPower.A3.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
% %             mean(tPower.A4.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan');...
% %             mean(tPower.S1.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
% %             mean(tPower.S2.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
% %             mean(tPower.S3.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
% %             mean(tPower.S4.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan')];
% %         group_std = [std(tPower.A1.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
% %             std(tPower.A2.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
% %             std(tPower.A3.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
% %             std(tPower.A4.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan');...
% %             std(tPower.S1.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
% %             std(tPower.S2.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
% %             std(tPower.S3.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan'),...
% %             std(tPower.S4.meanfreqPC(:,band_i,chan_to_plot(chan_i)),'omitnan')];
% %         group_sem = group_std/sqrt(Num_sub_trimmed);
% % 
% %            
% %         subplot(ceil(sqrt(length(chan_to_plot))),ceil(sqrt(length(chan_to_plot))),chan_i);%(ceil(sqrt(length(chan_to_plot))),ceil(sqrt(length(chan_to_plot))),chan_i);%(6,10,chan_i)
% %         barwitherr(group_sem,group_mean)
% %         colormap(gray)
% %         set(gca,'XTickLabel',{'Active','Sham'})
% %         ylabel(['Mean Power ' Band_labels{band_i} ' (microV2) % Baseline'])
% %         title([Band_labels{band_i} ' in ' chan_to_plot_labels{chan_i}]);
% %         
% %     end % chan_i
% %     legend({'PretACS';'After 5min';'After 10min';'After 15min'})
% % end % band_i
% % 
%  
% %% Topographical maps of % baseline values
%  
% Plot the power values at baseline topographically per frequency band
% for band_i = 1:Num_band
%     limMaxA = max(tPower.A.meanmeanfreq(:,band_i,:));
%     limMaxS = max(tPower.S.meanmeanfreq(:,band_i,:));
%     limMax = max(limMaxA,limMaxS);
%     
%     figure
%     subplot(2,2,1)
%     topoplot(tPower.A1.meanmeanfreq(1,band_i,:),EEG.chanlocs,'maplimits',[0 limMax],'electrodes','labels','conv','on','plotchans', chan_to_plot);
%     cbar('vert',0,[0 limMax]);
%     %cbar('vert',0,[-1 1]*max(abs(tPower.A.meanmeanfreqPC(1,band_i,:))));
%     %legend({'Active'})
%     title(['Meanfreq power (uV2) at baseline in ' Band_labels(band_i)]);
%     subplot(2,2,2)
%     topoplot(tPower.S1.meanmeanfreq(1,band_i,:),EEG.chanlocs,'maplimits',[0 limMax],'electrodes','labels','conv','on','plotchans', chan_to_plot);
%     cbar('vert',0,[0 limMax]);
%     %cbar('vert',0,[-1 1]*max(abs(tPower.S.meanmeanfreqPC(1,band_i,:))));
%     %legend({'Sham'})
%     subplot(2,2,3)
%     topoplot(tPower.A.meanmeanfreq(1,band_i,:),EEG.chanlocs,'maplimits',[0 limMax],'electrodes','labels','conv','on','plotchans', chan_to_plot);
%     cbar('vert',0,[0 limMax]);
%     subplot(2,2,4)
%     topoplot(tPower.S.meanmeanfreq(1,band_i,:),EEG.chanlocs,'maplimits',[0 limMax],'electrodes','labels','conv','on','plotchans', chan_to_plot);
%     cbar('vert',0,[0 limMax]);
% end %band_i
% 
% % Plot the power values topographically per frequency band
% for band_i = 1:Num_band
%     limMaxA = max(tPower.A.meanmeanfreq(:,band_i,:));
%     limMaxS = max(tPower.S.meanmeanfreq(:,band_i,:));
%     limMax = max(limMaxA,limMaxS);
%     
%     figure
%     subplot(1,2,1)
%     topoplot(tPower.A1.meanmeanfreq(1,band_i,:),EEG.chanlocs,'maplimits',[0 limMax],'electrodes','labels','conv','on','plotchans', chan_to_plot);
%     cbar('vert',0,[0 limMax]);
%     %cbar('vert',0,[-1 1]*max(abs(tPower.A.meanmeanfreqPC(1,band_i,:))));
%     %legend({'Active'})
%     title(['Meanfreq power (uV2) in ' Band_labels(band_i)]);
%     subplot(1,2,2)
%     topoplot(tPower.S1.meanmeanfreq(1,band_i,:),EEG.chanlocs,'maplimits',[0 limMax],'electrodes','labels','conv','on','plotchans', chan_to_plot);
%     cbar('vert',0,[0 limMax]);
%     %cbar('vert',0,[-1 1]*max(abs(tPower.S.meanmeanfreqPC(1,band_i,:))));
%     %legend({'Sham'})
%     
% end %band_i

% 
% % % Plot the % change topographically per subject (won't work for NaNs...)
% for sub_i = 1:Num_sub
% %sub_i= 2; % whomever you want to see!
% % Knowing that Band_values = [0.50,1;1,4;4,8;8,15;15,30;30,60]; % min and max (columns) for each frequency band of interest (rows)
% % and Band_labels = {'tACS frequency', 'delta', 'theta', 'alpha', 'beta', 'gamma'};
% 
% if isnan(tPower.A.meanfreq(sub_i,1,1))
% disp('NaN')
%     
% else limMaxA = max(abs(tPower.A.meanfreq_dB(sub_i,1,:)));
% limMaxS = max(abs(tPower.S.meanfreq_dB(sub_i,1,:)));
% limMax = max(limMaxA,limMaxS);
% 
%     figure
%     subplot(1,2,1)
%     topoplot(tPower.A.meanfreq_dB(sub_i,1,:),EEG.chanlocs,'electrodes','labels')
%     cbar('vert',0,[0 limMax]);
%     legend({'Active'})
%     title(['MeanfreqdB in tACS freq for subject S' num2str(sub_i)]);
%     subplot(1,2,2)
%     topoplot(tPower.S.meanfreq_dB(sub_i,1,:),EEG.chanlocs,'electrodes','labels')
%     cbar('vert',0,[0 limMax]);
%     legend({'Sham'})
% end
% end %sub_i

