% power analysis
% plot demo: power in different frequency band
% lzr
clear; clc; 
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StimSite = 'L_DLPFC';
RecSite = 'amygdala';
StimFreq = '10Hz';
Sham_cutaneous = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['rTMS_power_L_DLPFC_10Hz_' RecSite '.mat']);

if strcmp(StimSite, 'L_DLPFC')
    if strcmp(StimFreq, '10Hz')
        if strcmp(RecSite, 'amygdala')
            if Sham_cutaneous == 0
                channel_list = {'430_LFPx43' '430_LFPx44' '430_LFPx8' '430_LFPx9' '430_LFPx10' '430_LFPx11' '430_LFPx12' '430_LFPx13' ...
                                '460_LFPx105' '460_LFPx106' '460_LFPx107' ...
                                '524_LFPx149' ...
                                '534_LFPx231' '534_LFPx229' '534_LFPx230' '534_LFPx228'};
            end
        end
    end
end

%%
% from your power analysis: 
% the time point and freq information for the matrix
times = Power_times;
freqs = Power_freqs;

mask1 = ~(freqs > 55 & freqs <= 65); % Mask for values outside [55, 65]
mask2 = ~(freqs > 115 & freqs <= 125); % Mask for values outside [115, 125]
final_mask = mask1 & mask2;

freqs = freqs(final_mask);


% original data
% organized: channel*freqs*time
tms_data = zeros(height(Power_channel_tms), length(Power_freqs), length(Power_times));
for i = 1:height(Power_channel_tms)
    tms_data(i, :, :) = Power_channel_tms{i};
end

sham_data = zeros(height(Power_channel_sham), length(Power_freqs), length(Power_times));
for i = 1:height(Power_channel_sham)
    sham_data(i, :, :) = Power_channel_sham{i};
end

fs = 1000;

% the time before stimulation onset
before_time = 5; 

% the time range for plot
draw_time = [-2,9];
draw_pos = (draw_time)*fs;

% the time range for fdr
% attention: as "times" was not corrected as onset, the fdr_time here should
% also match the time of "times"
fdr_time = [4,9];
fdr_pos = (fdr_time + before_time)*fs;


% frequency bands
frequency_bands = {};
frequency_bands{1,1} = [4,8];
frequency_bands{2,1} = [8,12];
frequency_bands{3,1} = [13,30];
frequency_bands{4,1} = [30,50];
frequency_bands{5,1} = [70,110];
frequency_name = {};
frequency_name{1,1} = 'Theta';
frequency_name{2,1} = 'Alpha';
frequency_name{3,1} = 'Beta';
frequency_name{4,1} = 'Gamma';
frequency_name{5,1} = 'High-Gamma';
fre_num = size(frequency_bands,1);

frequency_bands = {};
frequency_bands{1,1} = [4,8];
frequency_bands{2,1} = [8,12];
frequency_bands{3,1} = [13,30];
frequency_bands{4,1} = [30,150];
frequency_name = {};
frequency_name{1,1} = 'Theta';
frequency_name{2,1} = 'Alpha';
frequency_name{3,1} = 'Beta';
frequency_name{4,1} = 'Broadband Gamma';
fre_num = size(frequency_bands,1);

sig_pos_freqband = cell(1,length(frequency_name));

%%

figure;
set(gcf,'color',[1,1,1])
set(gcf,'outerposition',[100,100,2000,2000]);



for fre_i = 1:fre_num
    
    fre1 = frequency_bands{fre_i,1}(1);
    fre2 = frequency_bands{fre_i,1}(2);
    fre_pos = find(freqs >= fre1 & freqs < fre2);
    
    % ttest
    target_data = squeeze(mean(tms_data(:,fre_pos,:),2));
    control_data = squeeze(mean(sham_data(:,fre_pos,:),2));
    temp_p_data = [];
    for time_i = 1:size(target_data,2)
        [h,p,ci,stats] = ttest(target_data(:,time_i)-control_data(:,time_i),0);
        temp_p_data(time_i,1) = p;
    end
    
    % fdr correction
    fdr_select_pos = find(times > fdr_pos(1) & times < fdr_pos(2));
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(temp_p_data(fdr_select_pos));
    temp_q_data = temp_p_data;
    temp_q_data(fdr_select_pos) = adj_p;
    
    
    % average the channels
    plot_data_ave = mean(target_data-control_data,1);
    plot_data_std = std(target_data-control_data,1)/sqrt(height(Power_channel_tms));
    
    
    subplot(fre_num,1,fre_i)
    time = times - before_time*fs;
    hold on
    plot(time, zeros(length(plot_data_ave),1),'Color',[1,1,1]*0.2,'linewidth',2)
    hold on
    ck_shadedErrorBar(time,plot_data_ave,plot_data_std,{'-','Color','r'},2);
    hold on
    
    % uncorrected significant
    hold on
    sig_pos = find(temp_p_data < 0.05);
    sig_pos = sig_pos(sig_pos > fdr_select_pos(1));
    plot(time(sig_pos),1*8,'*','Color',[1,1,1]*0.7,'linewidth',2)
    
    % fdr significant
    hold on
    sig_pos = find(temp_q_data < 0.05);
    sig_pos = sig_pos(sig_pos > fdr_select_pos(1));
    if ~isempty(sig_pos) 
        plot(time(sig_pos),0.93*8,'*','Color','r','linewidth',2)
    end
    sig_pos_freqband{fre_i} = sig_pos;
    
    % plot the zeros line
    hold on
    plot([0,0],[-8,8],'k:','linewidth',2)
    
    set(gca,'linewidth',2)
    set(gca,'Fontsize',15)
    
    xlim(draw_pos)
%     ylim([-1,1]*3)
    title(['TMS vs Sham, ' num2str(fre1),'-',num2str(fre2),'Hz'])
    
end

%%
power_sig_mean_table = {'Channel','theta_sig_tms','theta_sig_tvs','alpha_sig_tms','alpha_sig_tvs' ...
    'bgamma_sig_tms', 'bgamma_sig_tvs'};  
theta_tms_mean = cell(1,height(Power_channel_tms));
theta_tvs_mean = cell(1,height(Power_channel_tms));
alpha_tms_mean = cell(1,height(Power_channel_tms));
alpha_tvs_mean = cell(1,height(Power_channel_tms));
broadgamma_tms_mean = cell(1,height(Power_channel_tms));
broadgamma_tvs_mean = cell(1,height(Power_channel_tms));

rowIndex = 2;
for chl = 1:length(channel_list)
    % theta mean
    theta_tms_chl = mean(Power_channel_tms{chl,1}(:,sig_pos_freqband{1}));
    theta_sham_chl = mean(Power_channel_sham{chl,1}(:,sig_pos_freqband{1}));
    theta_tms_mean{chl} = mean(theta_tms_chl);
    theta_tvs_mean{chl} = mean(theta_tms_chl) - mean(theta_sham_chl);
    
    % alpha mean
    alpha_tms_chl = mean(Power_channel_tms{chl,1}(:,sig_pos_freqband{2}));
    alpha_sham_chl = mean(Power_channel_sham{chl,1}(:,sig_pos_freqband{2}));
    alpha_tms_mean{chl} = mean(alpha_tms_chl);
    alpha_tvs_mean{chl} = mean(alpha_tms_chl) - mean(alpha_sham_chl);
    
    % broadband gamma mean
    broadgamma_tms_chl = mean(Power_channel_tms{chl,1}(:,sig_pos_freqband{4}));
    broadgamma_sham_chl = mean(Power_channel_sham{chl,1}(:,sig_pos_freqband{4}));
    broadgamma_tms_mean{chl} = mean(broadgamma_tms_chl);
    broadgamma_tvs_mean{chl} = mean(broadgamma_tms_chl) - mean(broadgamma_sham_chl);
    
    power_sig_mean_table{rowIndex,1} = channel_list{chl};
    power_sig_mean_table{rowIndex,2} = theta_tms_mean{chl};
    power_sig_mean_table{rowIndex,3} = theta_tvs_mean{chl};
    power_sig_mean_table{rowIndex,4} = alpha_tms_mean{chl};
    power_sig_mean_table{rowIndex,5} = alpha_tvs_mean{chl};
    power_sig_mean_table{rowIndex,6} = broadgamma_tms_mean{chl};
    power_sig_mean_table{rowIndex,7} = broadgamma_tvs_mean{chl};
    rowIndex = rowIndex + 1;
end



writecell(power_sig_mean_table, ...
    '/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/3_rTMS_iEEG/4_DataAnalysis/Power_sig_mean.csv');











