%% TIME-FREQ ANALYSIS - for rTMS-iEEG
% Xianqing Bella Liu
% 12/04/2023
% Toolbox: EEGLAB

clear; clc; 
close all;

%%
rootpath = '/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/3_rTMS_iEEG';
savepath = [rootpath '/3_ProcessedData'];
cd(rootpath);

addpath /Users/xianqliu/Documents/MATLAB/fieldtrip-20230206
addpath /Users/xianqliu/Documents/MATLAB/eeglab2023.1
addpath(genpath('./0_Scripts'), ...
        genpath('./3_ProcessedData'));
    
ft_defaults;
eeglab;

%% 0. Initializing parameters and datasets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StimSite = 'L_DLPFC';
RecSite = 'hippocampus';
StimFreq = '10Hz';
Sham_cutaneous = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

savepath = [rootpath '/3_ProcessedData/' StimSite '_' StimFreq '_' RecSite];
mkdir(savepath);

sampling_rate = 1000;
preStimSecs = -5; % pre stim time in seconds (negative if before)
postStimSecs = 9; % post stim time in seconds
nTimepoint = (postStimSecs - preStimSecs)*sampling_rate + 1;

if strcmp(StimSite, 'L_DLPFC')
    if strcmp(StimFreq, '10Hz')
        if strcmp(RecSite, 'amygdala')
            if Sham_cutaneous == 0
                patient_list = [430 460 524 534];
                chlinfo = readtable([rootpath '/2_Info/Selected_channel.xlsx'], 'Sheet', 'rTMS L_DLPFC amyg');
                channel_list = cell(1,length(patient_list));
                channel_list{1} = {'LFPx43' 'LFPx44' 'LFPx8' 'LFPx9' 'LFPx10' 'LFPx11' 'LFPx12' 'LFPx13'}; % 430
                channel_list{2} = {'LFPx105' 'LFPx106' 'LFPx107'}; % 460
                channel_list{3} = {'LFPx149'}; % 524
                channel_list{4} = {'LFPx231' 'LFPx229' 'LFPx230' 'LFPx228'}; % 534
                nChannel = length([channel_list{:}]);
            end
        elseif strcmp(RecSite, 'hippocampus')
            if Sham_cutaneous == 0
                patient_list = [430 460 524];
                chlinfo = readtable([rootpath '/2_Info/Selected_channel.xlsx'], 'Sheet', 'rTMS L_DLPFC hipp');
                channel_list = cell(1,length(patient_list));
                channel_list{1} = {'LFPx51' 'LFPx52'}; % 430
                channel_list{2} = {'LFPx83' 'LFPx84' 'LFPx85' 'LFPx86' 'LFPx102' 'LFPx103'} ; % 460
                channel_list{3} = {'LFPx129' 'LFPx150' 'LFPx151' 'LFPx152' 'LFPx153' 'LFPx154' 'LFPx236'}; % 524
                nChannel = length([channel_list{:}]);
            end
        end
    end
end

%% 1. Load Data for Each Patient
% Create a cell array to hold the data from each patient
all_patients_data = cell(1, length(patient_list)); 
all_patients_tms = cell(1, length(patient_list));
all_patients_sham = cell(1, length(patient_list));

% Loop through each patient's dataset
for i = 1:length(patient_list)
    data = load([rootpath '/3_ProcessedData/' num2str(patient_list(i)) ...
        '_Comparison_' StimSite '_' StimFreq '_4A-CleanEpoched.mat']);
    
    all_patients_data{i} = data;
end

clear data ftData trigLengths trigShift trigTimes sessionDir


%% 2. Select Channels
% TMS
for i = 1:length(patient_list)
    cfg = [];
    cfg.channel = channel_list{i}; % amyg channel each patient

    all_patients_tms{i} = ft_selectdata(cfg, all_patients_data{i}.ftData_epoch_tms);
end

% Sham
for i = 1:length(patient_list)
    cfg = [];
    cfg.channel = channel_list{i}; % amyg channel each patient

    all_patients_sham{i} = ft_selectdata(cfg, all_patients_data{i}.ftData_epoch_sham);
end


%% 3. Select Trials
% TMS
all_patients_tms{1,3}.trial(:, 4) = []; % 524
all_patients_tms{1,3}.time(:, 4) = []; 

% Sham
all_patients_sham{1,3}.trial(:, 2) = []; % 524
all_patients_sham{1,3}.time(:, 2) = []; 


%% 4. Prepare Data for power analysis
% all_patients_sham{1, i_patient}.trial{1, i_trial}(i_channel, i_time) -->
% all_sham{1, i_patient}(i_channel, i_time, i_trial)

% TMS
all_tms = cell(1, length(patient_list));
for i_patient = 1:length(patient_list)
    all_tms{1,i_patient} = zeros(length(channel_list{i_patient}), nTimepoint, length(all_patients_tms{1, i_patient}.trial));
    for i_trial = 1:length(all_patients_tms{1, i_patient}.trial)
        for i_channel = 1:length(channel_list{i_patient})
            for i_time = 1:nTimepoint
                all_tms{1, i_patient}(i_channel, i_time, i_trial) = all_patients_tms{1, i_patient}.trial{1, i_trial}(i_channel, i_time);
            end
        end
    end
end

% Sham
all_sham = cell(1, length(patient_list));
for i_patient = 1:length(patient_list)
    all_sham{1,i_patient} = zeros(length(channel_list{i_patient}), nTimepoint, length(all_patients_sham{1, i_patient}.trial));
    for i_trial = 1:length(all_patients_sham{1, i_patient}.trial)
        for i_channel = 1:length(channel_list{i_patient})
            for i_time = 1:nTimepoint
                all_sham{1, i_patient}(i_channel, i_time, i_trial) = all_patients_sham{1, i_patient}.trial{1, i_trial}(i_channel, i_time);
            end
        end
    end
end

%% 5. Time-Frequency analysis

% set tf parameters
MarkTimes = [0 4000]; % e.g. [-3.25 -3 0 5]*1000; % this would plot a marker at the x zero crossing. add any other times in secs to also plot these
baseline = [-2000 -500];
baselineToRemove = [-2000 -500]+5000; % % baselineTime = [startTimeInMs endTimeInMs]; % in ms; removes this baseline prior to wavelet analysis
timesout = 2500; % temporal resolution of wavelet - divide entire duration (secs) of epoch by this to get resolution in seconds
alpha_val = 0.05; % significance threshold
padratio = 4;
waveFreq = [4 200]; % frequency range

% plotting details:
maxersp = 4;
itcplot = 'off';
ersplim = [-maxersp maxersp];
trialBaseStatus = 'on'; % either 'on', 'off' or 'full'

cycle_num = [2,0.5];

tf_parameter = struct('MarkTimes',MarkTimes,'baselineToRemove',baselineToRemove,'timesout',timesout,...
    'alpha_val',alpha_val,'preStimSecs',preStimSecs,'postStimSecs',postStimSecs,'padratio',padratio,...
    'waveFreq',waveFreq,'maxersp',maxersp,'itcplot',itcplot,'ersplim',ersplim,...
    'trialBaseStatus',trialBaseStatus,'cycle_num',cycle_num);
channel_tf_result = [];

% TMS
for i_patient = 1:length(patient_list)
    for i_channel = 1:length(channel_list{i_patient})
        tic
        data = [];
        data(1,:,:) = all_tms{1,i_patient}(i_channel,:,:); % 3-dimension, channel * time * trial
        
        NumTrials = size(data,3);
        
        
        % load data into EEG lab and clean
        EEG = pop_importdata( 'dataformat', 'array', 'data', 'data', 'setname', 'AllEpochs', 'srate',sampling_rate, 'pnts',0, 'nbchan',0);
        EEG = eeg_checkset( EEG );
        EEG = pop_importepoch(EEG, ones(1,NumTrials),{'AllEpochs'}, 'latencyfields',{'AllEpochs'}, 'timeunit',1, 'headerlines',0);
        EEG = eeg_checkset( EEG );
        
        
        currData = EEG;
        elecidx = 1;
        
        [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = pop_newtimef(currData, ...
            1, elecidx, [EEG.xmin EEG.xmax]*1000, cycle_num,'freqs',[waveFreq(1) waveFreq(2)],'padratio', padratio, ...
            'plotphase', 'off', 'timesout',timesout,'naccu', 200, 'alpha',alpha_val,'baseboot',1,'rmerp','off','freqscale','log', ...
            'plotersp','off', 'plotitc',itcplot,'marktimes',MarkTimes,'erspmax', maxersp,'trialbase',trialBaseStatus,'baseline',baselineToRemove);
        
        fddata_result2 = abs(tfdata);
        
%         channel_tf_result.ave = ersp;  % average tf among trials, it could be baseline correction depending on the 'trialBaseStatus'
        channel_tf_result.times = times;  % [times]: the returned time value for tf results
        channel_tf_result.freqs = freqs;  % [freqs]: the returned frequency value  for tf results
        channel_tf_result.itc = itc;   % PLV
        channel_tf_result.tfData = tfdata; % tf before trial average, complex value
        
%         sub = all_channel_data_raw.info{i_channel,1};
%         channel_info = all_channel_data_raw.info{i_channel,2};

        savename = [savepath '/' num2str(patient_list(i_patient)) '_' channel_list{i_patient}{i_channel} ...
            '_power_' StimSite '_' StimFreq '_' RecSite '_TMS.mat'];
        mkdir([savepath '/' StimSite '_' StimFreq '_' RecSite]);
        save(savename,'channel_tf_result','tf_parameter','-v7.3');
        
        temp_time = toc;
        %     disp(['------Working ', num2str(sub),'-',num2str(channel_info),':',num2str(temp_time),'s'])
    end
end

% Sham
for i_patient = 1:length(patient_list)
    for i_channel = 1:length(channel_list{i_patient})
        tic
        data = [];
        data(1,:,:) = all_sham{1,i_patient}(i_channel,:,:); % 3-dimension, channel * time * trial
        
        NumTrials = size(data,3);
        
        
        % load data into EEG lab and clean
        EEG = pop_importdata( 'dataformat', 'array', 'data', 'data', 'setname', 'AllEpochs', 'srate',sampling_rate, 'pnts',0, 'nbchan',0);
        EEG = eeg_checkset( EEG );
        EEG = pop_importepoch(EEG, ones(1,NumTrials),{'AllEpochs'}, 'latencyfields',{'AllEpochs'}, 'timeunit',1, 'headerlines',0);
        EEG = eeg_checkset( EEG );
        
        
        currData = EEG;
        elecidx = 1;
        
        [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = pop_newtimef(currData, ...
            1, elecidx, [EEG.xmin EEG.xmax]*1000, cycle_num,'freqs',[waveFreq(1) waveFreq(2)],'padratio', padratio, ...
            'plotphase', 'off', 'timesout',timesout,'naccu', 200, 'alpha',alpha_val,'baseboot',1,'rmerp','off','freqscale','log', ...
            'plotersp','off', 'plotitc',itcplot,'marktimes',MarkTimes,'erspmax', maxersp,'trialbase',trialBaseStatus,'baseline',baselineToRemove);
        
        fddata_result2 = abs(tfdata);
        
%         channel_tf_result.ave = ersp;  % average tf among trials, it could be baseline correction depending on the 'trialBaseStatus'
        channel_tf_result.times = times;  % [times]: the returned time value for tf results
        channel_tf_result.freqs = freqs;  % [freqs]: the returned frequency value  for tf results
        channel_tf_result.itc = itc;   % PLV
        channel_tf_result.tfData = tfdata; % tf before trial average, complex value
        
%         sub = all_channel_data_raw.info{i_channel,1};
%         channel_info = all_channel_data_raw.info{i_channel,2};
        
        savename = [savepath '/' num2str(patient_list(i_patient)) '_' channel_list{i_patient}{i_channel} ...
            '_power_' StimSite '_' StimFreq '_' RecSite '_Sham.mat'];
        save(savename,'channel_tf_result','tf_parameter','-v7.3');
        
        temp_time = toc;
        %     disp(['------Working ', num2str(sub),'-',num2str(channel_info),':',num2str(temp_time),'s'])
    end
end


%% 6. plot time frequency power spectrum
% savepath = [savepath '/Power'];
% mkdir(savepath);

% tms
files_tms = dir([savepath '/*TMS.mat']);
Power_channel_tms = cell(nChannel, 1);
for i_channel = 1:height(files_tms)
    
    load([savepath '/' files_tms(i_channel).name]);
    fprintf(['Working on ' files_tms(i_channel).name ' \n']);
    clear P

    alltfX = channel_tf_result.tfData;
    P  = alltfX.*conj(alltfX);
    P = 10 * log10(P);  % log scale
    
%     % baseline subtraction
%     baselineToRemove = [-2500 -500]+5000;
%     base_pos = find(channel_tf_result.times > min(baselineToRemove) & channel_tf_result.times < max(baselineToRemove));
%     P = P - mean(P(:,base_pos),2);
    
%     figure;
%     set(gcf,'color',[1,1,1])
%     set(gcf,'outerposition',[0,0,3000,1000]);
% 
%     for i_trial = 1:size(P, 3)
%         subplot(4,2,i_trial)
%         do_EEGLABPlot(P(:,:,i_trial), ...
%             [],channel_tf_result.times-5000,channel_tf_result.freqs,[-1,1]*12,[0],[files_tms(i_channel).name(1:end-4)],0)
%         set(gca,'fontsize',18)
% %         xlim(time_limits)
%     end
%     
%     savename = [savepath, '/' files_tms(i_channel).name(1:end-4),'.png'];
%     img = gcf;
%     print(img,'-dpng','-r450',savename)
    
    P_avg = mean(P, 3);
    
    % baseline subtraction
    base_pos = find(channel_tf_result.times > min(baselineToRemove) & channel_tf_result.times < max(baselineToRemove));
    P_avg = P_avg - mean(P_avg(:,base_pos),2);
    
    
%     figure;
%     set(gcf,'color',[1,1,1])
%     set(gcf,'outerposition',[0,0,2700,1000]);
%     
%     do_EEGLABPlot(P_avg(:,:,1), ...
%             [],channel_tf_result.times-5000,channel_tf_result.freqs,[-1,1]*12,[0],[files_tms(i_channel).name(1:end-4)],0)
%         set(gca,'fontsize',18)
%         
%     savename = [savepath, '/Power/' files_tms(i_channel).name(1:end-4),'_avg.png'];
%     img = gcf;
%     print(img,'-dpng','-r450',savename)
    
    Power_channel_tms{i_channel,1} = P_avg;
end

% sham
files_sham = dir([savepath '/*Sham.mat']);
Power_channel_sham = cell(nChannel, 1);
for i_channel = 1:height(files_sham)
    
    load([savepath '/' files_sham(i_channel).name]);
    fprintf(['Working on ' files_sham(i_channel).name ' \n']);
    clear P

    alltfX = channel_tf_result.tfData;
    P  = alltfX.*conj(alltfX);
    P = 10 * log10(P);  % log scale
    
%     % baseline subtraction
%     base_pos = find(channel_tf_result.times > min(baselineToRemove) & channel_tf_result.times < max(baselineToRemove));
%     P = P - mean(P(:,base_pos),2);
    
%     figure;
%     set(gcf,'color',[1,1,1])
%     set(gcf,'outerposition',[0,0,3000,1000]);
% 
%     for i_trial = 1:size(P, 3)
%         subplot(4,2,i_trial)
%         do_EEGLABPlot(P(:,:,i_trial), ...
%             [],channel_tf_result.times-5000,channel_tf_result.freqs,[-1,1]*12,[0],[files_sham(i_channel).name(1:end-4)],0)
%         set(gca,'fontsize',18)
% %         xlim(time_limits)
%     end
%     
%     savename = [savepath, '/' files_sham(i_channel).name(1:end-4),'.png'];
%     img = gcf;
%     print(img,'-dpng','-r450',savename)
    
    P_avg = mean(P, 3);
    
    % baseline subtraction
    base_pos = find(channel_tf_result.times > min(baselineToRemove) & channel_tf_result.times < max(baselineToRemove));
    P_avg = P_avg - mean(P_avg(:,base_pos),2);
    
%     figure;
%     set(gcf,'color',[1,1,1])
%     set(gcf,'outerposition',[0,0,2700,1000]);
%     
%     do_EEGLABPlot(P_avg(:,:,1), ...
%             [],channel_tf_result.times-5000,channel_tf_result.freqs,[-1,1]*12,[0],[files_sham(i_channel).name(1:end-4)],0)
%         set(gca,'fontsize',18)
%         
%     savename = [savepath, '/Power/' files_sham(i_channel).name(1:end-4),'_avg.png'];
%     img = gcf;
%     print(img,'-dpng','-r450',savename)
    
    Power_channel_sham{i_channel,1} = P_avg;
end

Power_times = channel_tf_result.times;
Power_freqs = channel_tf_result.freqs;


%% TMS vs Sham t statistics

% Initialize t-value array
Power_t_values_tms = zeros(size(Power_channel_tms{1}));
Power_p_values_tms = zeros(size(Power_channel_tms{1}));
Power_t_values_sham = zeros(size(Power_channel_tms{1}));
Power_p_values_sham = zeros(size(Power_channel_tms{1}));
Power_t_values = zeros(size(Power_channel_tms{1}));
Power_p_values = zeros(size(Power_channel_tms{1}));

% Loop over each frequency and time point
for freq = 1:height(Power_channel_tms{1})
    for time = 1:width(Power_channel_tms{1})
        % Extract data for the current frequency and time across all channels
        tms_values = cellfun(@(x) x(freq, time), Power_channel_tms);
        sham_values = cellfun(@(x) x(freq, time), Power_channel_sham);
        
        % one-sample t test
        [h, p, ~, stats] = ttest(tms_values, 0);
        Power_t_values_tms(freq, time) = stats.tstat;
        Power_p_values_tms(freq, time) = p; % Store the p-value
        
        [h, p, ~, stats] = ttest(sham_values, 0);
        Power_t_values_sham(freq, time) = stats.tstat;
        Power_p_values_sham(freq, time) = p; 

        % Paired t-test
        [h, p, ~, stats] = ttest(tms_values, sham_values);
        Power_t_values(freq, time) = stats.tstat;
        Power_p_values(freq, time) = p; 
    end
    fprintf(['Done with freq ' num2str(freq) ' \n']);
end


save([savepath '/rTMS_power_' StimSite '_' StimFreq '_' RecSite '.mat'],'Power_channel_tms','Power_channel_sham',...
    'Power_times','Power_freqs','Power_t_values_tms','Power_p_values_tms','Power_t_values_sham','Power_p_values_sham', ...
    'Power_t_values','Power_p_values','-v7.3');


%%  plot TMS, Sham, and TMS vs Sham t map
% TMS
figure;
set(gcf,'color',[1,1,1])
set(gcf,'outerposition',[0,0,2000,1000]);

subplot(2,1,1)
do_EEGLABPlot(Power_t_values_tms, [], Power_times-5000, Power_freqs,[-1,1]*6,...
    [0], ['TMS t value'],0)
set(gca,'fontsize',18)

subplot(2,1,2)
do_EEGLABPlot_masked(Power_t_values_tms, [], Power_times-5000, Power_freqs,[-1,1]*6,...
    [0], Power_p_values_tms, ['TMS t value (thresholded)'],0)
set(gca,'fontsize',18)

rootpath = '/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/3_rTMS_iEEG';
figurepath = [rootpath '/1_Figures/Power'];
savename = [figurepath, '/TMS_' StimSite '_' StimFreq '_' RecSite '_t_sig.png'];
img = gcf;
print(img,'-dpng','-r450',savename)


% Sham
figure;
set(gcf,'color',[1,1,1])
set(gcf,'outerposition',[0,0,2000,1000]);

subplot(2,1,1)
do_EEGLABPlot(Power_t_values_sham, [], Power_times-5000, Power_freqs,[-1,1]*6,...
    [0], ['Sham t value'],0)
set(gca,'fontsize',18)

subplot(2,1,2)
do_EEGLABPlot_masked(Power_t_values_sham, [], Power_times-5000, Power_freqs,[-1,1]*6,...
    [0], Power_p_values_sham, ['Sham t value (thresholded)'],0)
set(gca,'fontsize',18)

rootpath = '/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/3_rTMS_iEEG';
figurepath = [rootpath '/1_Figures/Power'];
savename = [figurepath, '/Sham_' StimSite '_' StimFreq '_' RecSite '_t_sig.png'];
img = gcf;
print(img,'-dpng','-r450',savename)


% TMS vs Sham
figure;
set(gcf,'color',[1,1,1])
set(gcf,'outerposition',[0,0,2000,1000]);

subplot(2,1,1)
do_EEGLABPlot(Power_t_values, [], Power_times-5000, Power_freqs,[-1,1]*6,...
    [0], ['TMS vs Sham t value'],0)
set(gca,'fontsize',18)

subplot(2,1,2)
do_EEGLABPlot_masked(Power_t_values, [], Power_times-5000, Power_freqs,[-1,1]*6,...
    [0], Power_p_values, ['TMS vs Sham t value (thresholded)'],0)
set(gca,'fontsize',18)

rootpath = '/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/3_rTMS_iEEG';
figurepath = [rootpath '/1_Figures/Power'];
savename = [figurepath, '/TMS_vs_Sham_' StimSite '_' StimFreq '_' RecSite '_t_sig.png'];
img = gcf;
print(img,'-dpng','-r450',savename)
























