
% Toolbox: EEGLAB

clear; clc; 
close all;

rootpath = '/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/3_rTMS_iEEG';
savepath = [rootpath '/3_ProcessedData/Power'];
cd(rootpath);

addpath /Users/xianqliu/Documents/MATLAB/fieldtrip-20230206
addpath /Users/xianqliu/Documents/MATLAB/eeglab2023.1
addpath(genpath('./0_Scripts'), ...
        genpath('./3_ProcessedData'));
    
ft_defaults;
eeglab;

%% Setup parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StimSite = 'L_DLPFC';
StimFreq = '10Hz';
RecSite = 'hippocampus';
Sham_cutaneous = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampling_rate = 1000;
preStimSecs = -5; % pre stim time in seconds (negative if before)
postStimSecs = 9; % post stim time in seconds
nTimepoint = (postStimSecs - preStimSecs)*sampling_rate + 1;
baseline = [-2500 500];

if strcmp(StimSite, 'L_DLPFC')
    if strcmp(StimFreq, '10Hz')
        if Sham_cutaneous == 0
            patient_list = [430 460 524 534];
            chlinfo = readtable([rootpath '/2_Info/Selected_channel_AMY.xlsx'], 'Sheet', 'rTMS L_DLPFC');
            amygChannel{1} = {'LFPx43' 'LFPx44' 'LFPx8' 'LFPx9' 'LFPx10' 'LFPx11' 'LFPx12' 'LFPx13'}; % 430
            amygChannel{2} = {'LFPx105' 'LFPx106' 'LFPx107'}; % 460
            amygChannel{3} = {'LFPx149'}; % 524
            amygChannel{4} = {'LFPx231' 'LFPx229' 'LFPx230' 'LFPx228'}; % 534
            nChannel = length([amygChannel{:}]);
        end
    end
end

%% plot time frequency power spectrum

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
    
    % baseline subtraction
    base_pos = find(channel_tf_result.times > min(baseline) & channel_tf_result.times < max(baseline));
    P = P - mean(P(:,base_pos),2);
    
    figure;
    set(gcf,'color',[1,1,1])
    set(gcf,'outerposition',[0,0,3000,1000]);

    for i_trial = 1:size(P, 3)
        subplot(4,2,i_trial)
        do_EEGLABPlot(P(:,:,i_trial), ...
            [],channel_tf_result.times-5000,channel_tf_result.freqs,[-1,1]*12,[0],[files_tms(i_channel).name(1:end-4)],0)
        set(gca,'fontsize',18)
%         xlim(time_limits)
    end
    
    savename = [savepath, '/' files_tms(i_channel).name(1:end-4),'.png'];
    img = gcf;
    print(img,'-dpng','-r450',savename)
    
    P_avg = mean(P, 3);
    
    figure;
    set(gcf,'color',[1,1,1])
    set(gcf,'outerposition',[0,0,2700,1000]);
    
    do_EEGLABPlot(P_avg(:,:,1), ...
            [],channel_tf_result.times-5000,channel_tf_result.freqs,[-1,1]*12,[0],[files_tms(i_channel).name(1:end-4)],0)
        set(gca,'fontsize',18)
        
    savename = [savepath, '/' files_tms(i_channel).name(1:end-4),'_avg.png'];
    img = gcf;
    print(img,'-dpng','-r450',savename)
    
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
    
    % baseline subtraction
    base_pos = find(channel_tf_result.times > min(baseline) & channel_tf_result.times < max(baseline));
    P = P - mean(P(:,base_pos),2);
    
    figure;
    set(gcf,'color',[1,1,1])
    set(gcf,'outerposition',[0,0,3000,1000]);

    for i_trial = 1:size(P, 3)
        subplot(4,2,i_trial)
        do_EEGLABPlot(P(:,:,i_trial), ...
            [],channel_tf_result.times-5000,channel_tf_result.freqs,[-1,1]*12,[0],[files_sham(i_channel).name(1:end-4)],0)
        set(gca,'fontsize',18)
%         xlim(time_limits)
    end
    
    savename = [savepath, '/' files_sham(i_channel).name(1:end-4),'.png'];
    img = gcf;
    print(img,'-dpng','-r450',savename)
    
    P_avg = mean(P, 3);
    
    figure;
    set(gcf,'color',[1,1,1])
    set(gcf,'outerposition',[0,0,2700,1000]);
    
    do_EEGLABPlot(P_avg(:,:,1), ...
            [],channel_tf_result.times-5000,channel_tf_result.freqs,[-1,1]*12,[0],[files_sham(i_channel).name(1:end-4)],0)
        set(gca,'fontsize',18)
        
    savename = [savepath, '/' files_sham(i_channel).name(1:end-4),'_avg.png'];
    img = gcf;
    print(img,'-dpng','-r450',savename)
    
    Power_channel_sham{i_channel,1} = P_avg;
end

Power_times = channel_tf_result.times;
Power_freqs = channel_tf_result.freqs;


%% TMS vs Sham t statistics

% Initialize t-value array
Power_t_values = zeros(size(Power_channel_tms{1}));
Power_p_values = zeros(size(Power_channel_tms{1}));

% Loop over each frequency and time point
for freq = 1:height(Power_channel_tms{1})
    for time = 1:width(Power_channel_tms{1})
        % Extract data for the current frequency and time across all channels
        tms_values = cellfun(@(x) x(freq, time), Power_channel_tms);
        sham_values = cellfun(@(x) x(freq, time), Power_channel_sham);

        % Paired t-test
        [~, ~, ~, stats] = ttest(tms_values, sham_values);
        [h, p, ~, stats] = ttest(tms_values, sham_values);

        % Store the t-value
        Power_t_values(freq, time) = stats.tstat;
        Power_p_values(freq, time) = p; % Store the p-value
    end
    fprintf(['Done with freq ' num2str(freq) ' \n']);
end


save([savepath '/rTMS_power_L_DLPFC_10Hz.mat'],'Power_channel_tms','Power_channel_sham',...
    'Power_times','Power_freqs','Power_t_values','Power_p_values','-v7.3');


%%  plot TMS vs Sham t map
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
savename = [figurepath, '/TMS_vs_Sham_t_sig.png'];
img = gcf;
print(img,'-dpng','-r450',savename)


%%








