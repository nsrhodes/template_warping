% Plot peak from desync for template warping comparison
%% Housekeeping
clear all
clc
close all

addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
ft_defaults;

%% set directories and sub numbers
project_dir =  'R:\DRS-KidsOPM\Temp_warp_paper\pseudoMRI\';
project_dir_indiv = 'R:\DRS-KidsOPM\Temp_warp_paper\individual\';
datadir = [project_dir,'Data',filesep,'BIDS',filesep];
datadir_indiv = [project_dir_indiv,'Data',filesep,'BIDS',filesep];

extension_T = '_4mm';
good_subs = [1:13 15 16 18 21 24:26];
path_mni_mesh = 'R:\DRS-KidsOPM\';
files_mni_mesh = 'meshes_mni.mat';
load([path_mni_mesh files_mni_mesh]); 

mni_brain = ft_read_mri([datadir, 'derivatives',filesep,'sourcespace',filesep,'MNI152_T1_1mm_brain.nii.gz']);
brain_mask = find(mni_brain.anatomy>0);

aal_regs = ft_read_mri([datadir, 'derivatives',filesep,'sourcespace',filesep,'AAL_regions_1mm_mas.nii.gz']);


for i = 1:78
    aal_brain = aal_regs.anatomy(:,:,:,i);
    aal_brain = aal_brain(brain_mask);
    aal_inds{i} = find(aal_brain>0);
end 

%% Loop through subjects loading in tstats in MNI space
count = 0;
for sub_i = good_subs
    count = count+1;
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'

    path_Tstat = [datadir,'derivatives',filesep,'Tstats',extension_T,filesep,'sub-',sub,filesep];
    path_Tstat_indiv = [datadir_indiv,'derivatives',filesep,'Tstats',extension_T,filesep,'sub-',sub,filesep];

    D1_tstat = ft_read_mri([path_Tstat 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_index_MNI.nii.gz']);
    D1_tstat_indiv = ft_read_mri([path_Tstat_indiv 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_index_MNI.nii.gz']);
    D4_tstat = ft_read_mri([path_Tstat 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_pinky_MNI.nii.gz']);
    D4_tstat_indiv = ft_read_mri([path_Tstat_indiv 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_pinky_MNI.nii.gz']);
    
    All_D1_tstats(:,:,:,count) = D1_tstat.anatomy;%./abs(min(D1_tstat.anatomy(:)));
    All_D1_tstats_indiv(:,:,:,count) = D1_tstat_indiv.anatomy;%./abs(min(D1_tstat_indiv.anatomy(:)));
    All_D4_tstats(:,:,:,count) = D4_tstat.anatomy./abs(min(D4_tstat.anatomy(:)));
    All_D4_tstats_indiv(:,:,:,count) = D4_tstat_indiv.anatomy./abs(min(D4_tstat_indiv.anatomy(:)));
end

%%
All_D1_tstats_vect = reshape(All_D1_tstats,[],20);
All_D1_tstats_indiv_vect = reshape(All_D1_tstats_indiv,[],20);

All_D1_tstats_vect = All_D1_tstats_vect(brain_mask,:);
All_D1_tstats_indiv_vect = All_D1_tstats_indiv_vect(brain_mask,:);

All_D4_tstats_vect = reshape(All_D4_tstats,[],20);
All_D4_tstats_indiv_vect = reshape(All_D4_tstats_indiv,[],20);

All_D4_tstats_vect = All_D4_tstats_vect(brain_mask,:);
All_D4_tstats_indiv_vect = All_D4_tstats_indiv_vect(brain_mask,:);
%% Correlate tstats
for i = 1:20
    [r,p] = corrcoef(All_D1_tstats_indiv_vect(:,i),All_D1_tstats_vect(:,i));
    tstat_cor(i) = r(2,1);
    [r_d4,p_d4] = corrcoef(All_D4_tstats_indiv_vect(:,i),All_D4_tstats_vect(:,i));
    tstat_cor_D4(i) = r_d4(2,1);
end 

% voilin plot 
addpath('C:\Users\ppanr\Documents\GitHub\Braille_kids\Violinplot-Matlab-master')
figure
violins = violinplot(tstat_cor);
ylim([0.6 1]); ylabel('Correlation between individual and pseudo-MRI'); xlim([0.6 1.4])


for i = 1:78
    inds = aal_inds{i};
    aal_tstat(i,:) = mean(All_D1_tstats_vect(inds,:),1);
    aal_tstat_indiv(i,:) = mean(All_D1_tstats_indiv_vect(inds,:),1);
end

aal_mean_diff = mean(aal_tstat-aal_tstat_indiv,2);
figure;
plot(mean(aal_tstat_indiv,2),mean(aal_tstat,2),'ko');
axis equal
hold on;plot([-0.3:0.1:0.1],[-0.3:0.1:0.1],'k--')
xlabel('Individual MRI pseudo-T'),ylabel('Pseudo-MRI pseudo-T');
xlim([-0.3 0.1]);
ylim([-0.3 0.1])

[aal_tstat_cor,p] = corr(mean(aal_tstat_indiv,2),mean(aal_tstat,2))

addpath('C:\Users\ppanr\Documents\GitHub\Braille_kids\BrainPlots')
addpath('C:\Users\ppanr\Documents\GitHub\Braille_kids\gifti-1.8\')
figure; subplot(1,3,1)
PaintBrodmannAreas_chooseview(aal_mean_diff, 78, 78, [-0.1 0.1], 'Pseudo-T',[],[90 0 0])

subplot(1,3,2)
PaintBrodmannAreas_chooseview(aal_mean_diff, 78, 78, [-0.1 0.1], 'Pseudo-T',[],[0 90 0])
subplot(1,3,3)
PaintBrodmannAreas_chooseview(aal_mean_diff, 78, 78, [-0.1 0.1], 'Pseudo-T',[],[0 0 90])
%colorbar
%%
% Save out average tstat

means = (All_D1_tstats_vect + All_D1_tstats_indiv_vect)./2;
diffs = (All_D1_tstats_vect - All_D1_tstats_indiv_vect);

mean_tstat = D1_tstat;
mean_tstat.anatomy = mean(All_D1_tstats,4);
mean_tstat_indiv = D1_tstat_indiv;
mean_tstat_indiv.anatomy = mean(All_D1_tstats_indiv,4);
mean_tstat_D4 = D4_tstat;
mean_tstat_D4.anatomy = mean(All_D4_tstats,4);
mean_tstat_indiv_D4 = D4_tstat_indiv;
mean_tstat_indiv_D4.anatomy = mean(All_D4_tstats_indiv,4);
ft_write_mri('R:\DRS-KidsOPM\Temp_warp_paper\mean_tstat-pseudoMRI',mean_tstat,'dataformat','nifti')
ft_write_mri('R:\DRS-KidsOPM\Temp_warp_paper\mean_tstat_indiv',mean_tstat_indiv,'dataformat','nifti')
ft_write_mri('R:\DRS-KidsOPM\Temp_warp_paper\mean_tstat-pseudoMRI_little',mean_tstat_D4,'dataformat','nifti')
ft_write_mri('R:\DRS-KidsOPM\Temp_warp_paper\mean_tstat_indiv_little',mean_tstat_indiv_D4,'dataformat','nifti')

[ra,pa] = corr(mean_tstat_indiv.anatomy(brain_mask),mean_tstat.anatomy(brain_mask))
[ra_d4,pa_d4] = corr(mean_tstat_indiv_D4.anatomy(brain_mask),mean_tstat_D4.anatomy(brain_mask))
%% 
% Plot average TFS and get channel and trial counts
count = 0;
fs = 1200;
conwin = [2.5,3]*fs;


addpath('C:\Users\ppanr\Documents\GitHub\Braille_kids\')
for sub_i = good_subs
    count = count+1;
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
    path_VE = [datadir,'derivatives',filesep,'Tstats',filesep,'sub-',sub,filesep];
    files_VE = ['D1_peak_VE_1-150_Hz.mat'];
    load([path_VE files_VE],'VE_index');
    trial_time = linspace(0,size(VE_index,1)./fs,size(VE_index,1));
    [TFS(:,:,count),fre] = VE_TFS_tw(VE_index,conwin,trial_time,fs);
    good_trials(count) = size(VE_index,2);
    load(['R:\DRS-KidsOPM\Paediatric_OPM_Notts_AdultData\Data\BIDS\derivatives\ICA\sub-' sub '\sub-' sub '_ses-001_task-braille_run-001_ICA_data.mat'])
    no_ch(count)=size(comp150.topolabel,1);
end
TFS_mean = mean(TFS,3);
figure
pcolor(trial_time,fre,TFS_mean);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
colorbar;caxis([-0.5 0.5])
axis fill; axis square; box off;
% Calc mean and sd trials 
mean_trials_removed = 41 - mean(good_trials);
sd_trials_removed = std(41-good_trials);
% Calc mean and sd channels
mean_ch_removed = 174-mean(no_ch);
sd_ch_removed = std(174-no_ch);
%% Plot average TFS individual and get channel and trial counts
count = 0;

addpath('C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\Braille_kids\')
for sub_i = good_subs
    count = count+1;
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
    path_VE = [datadir_indiv,'derivatives',filesep,'Tstats',filesep,'sub-',sub,filesep];
    files_VE = ['D1_peak_VE_1-150_Hz.mat'];
    load([path_VE files_VE],'VE_index');
    trial_time = linspace(0,size(VE_index,1)./fs,size(VE_index,1));
    [TFS_indiv(:,:,count),fre] = VE_TFS_tw(VE_index,conwin,trial_time,fs);
end
TFS_mean_indiv = mean(TFS_indiv,3);
figure
pcolor(trial_time,fre,TFS_mean_indiv);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
colorbar;caxis([-0.5 0.5])
axis fill; axis square; box off;

%% TFS correlation
for i = 1:20
    this_TFS = TFS(2:26,:,i);
    this_TFS_indiv = TFS_indiv(2:26,:,i);
    [r,p] = corrcoef(this_TFS(:),this_TFS_indiv(:));
    TFS_cor(i) = r(2,1);
end 

% voilin plot 
addpath('C:\Users\ppynr2\OneDrive - The University of Nottingham\matlab\tools\Violinplot-Matlab-master')
figure
violins = violinplot(TFS_cor');
ylim([0.5 1]); ylabel('Correlation between individual and pseudo-MRI'); xlim([0.6 1.4])
xticks('');

%%plot single TFS
figure
pcolor(trial_time,fre,TFS_indiv(:,:,1));shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
colorbar;caxis([-0.5 0.5])
axis fill; axis square; box off;

%% mean TFS correlation
TFS_mean(1,:) = [];
TFS_mean_indiv(1,:) = [];
[TFS_mean_cor,~] = corr(TFS_mean(:),TFS_mean_indiv(:))

%% D4
count = 0
for sub_i = good_subs
    count = count+1;
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
    path_VE = [datadir,'derivatives',filesep,'Tstats',filesep,'sub-',sub,filesep];
    files_VE = ['D4_peak_VE_1-150_Hz.mat'];
    load([path_VE files_VE],'VE_pinky');
    [TFS_D4(:,:,count),fre] = VE_TFS_tw(VE_pinky,conwin,trial_time,fs);
end
TFS_mean_D4 = mean(TFS_D4,3);
figure
pcolor(trial_time,fre,TFS_mean_D4);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
colorbar;caxis([-0.5 0.5])
axis fill; axis square; box off;

%% Plot average TFS individual and get channel and trial counts
count = 0;

addpath('C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\Braille_kids\')
for sub_i = good_subs
    count = count+1;
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
    path_VE = [datadir_indiv,'derivatives',filesep,'Tstats',filesep,'sub-',sub,filesep];
    files_VE = ['D4_peak_VE_1-150_Hz.mat'];
    load([path_VE files_VE],'VE_pinky');
    [TFS_indiv_D4(:,:,count),fre] = VE_TFS_tw(VE_pinky,conwin,trial_time,fs);
end
TFS_mean_indiv_D4 = mean(TFS_indiv_D4,3);
figure
pcolor(trial_time,fre,TFS_mean_indiv_D4);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
colorbar;caxis([-0.5 0.5])
axis fill; axis square; box off;

%% TFS correlation
for i = 1:20
    this_TFS = TFS_D4(2:26,:,i);
    this_TFS_indiv = TFS_indiv_D4(2:26,:,i);
    [r,p] = corrcoef(this_TFS(:),this_TFS_indiv(:));
    TFS_cor_D4(i) = r(2,1);
end 

%% mean TFS correlation
TFS_mean_D4(1,:) = [];
TFS_mean_indiv_D4(1,:) = [];
[TFS_mean_cor_D4,~] = corr(TFS_mean_D4(:),TFS_mean_indiv_D4(:))

%% plot single TFS
select = 3; 

figure
subplot(2,2,1)
pcolor(trial_time,fre,TFS(:,:,select));shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
colorbar;caxis([-0.5 0.5])
axis fill; axis square; box off;
set(gca,'FontSize',12,'FontName','Arial')
%title('pseudoMRI index')

subplot(2,2,2)
pcolor(trial_time,fre,TFS_indiv(:,:,select));shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
colorbar;caxis([-0.5 0.5])
axis fill; axis square; box off;
%title('individual index')
set(gca,'FontSize',12,'FontName','Arial')

subplot(2,2,3)
pcolor(trial_time,fre,TFS_D4(:,:,select));shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
colorbar;caxis([-0.5 0.5])
axis fill; axis square; box off;
%title('pseudoMRI little')
set(gca,'FontSize',12,'FontName','Arial')

subplot(2,2,4)
pcolor(trial_time,fre,TFS_indiv_D4(:,:,select));shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
colorbar;caxis([-0.5 0.5])
axis fill; axis square; box off;
%title('individual little')
set(gca,'FontSize',12,'FontName','Arial')

%% Gaussian over hist

% Generate 20 random data points from a normal distribution
data = tstat_cor;  % You can modify this line to generate different data

% Create a histogram
figure;
histogram(data,3, 'Normalization', 'pdf');  % Normalize to get probability density

% Hold the current figure to overlay the Gaussian fit
hold on;

% Fit a Gaussian to the data
mu = mean(data);        % Mean of the data
sigma = std(data);      % Standard deviation of the data

% Create a range of x values for the Gaussian curve
x = linspace(min(data)-1, max(data)+1, 100);
gaussianFit = (1/(sigma * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sigma).^2);

% Plot the Gaussian fit
plot(x, gaussianFit, 'r-', 'LineWidth', 2);

% Add titles and labels
title('Histogram with Gaussian Fit');
xlabel('Data Points');
ylabel('Probability Density');
legend('Histogram', 'Gaussian Fit');

% Release the hold
hold off;
