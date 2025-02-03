% Plot peak from desync for template warping comparison
%% Housekeeping
clear all
clc
close all

addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
addpath('C:\Users\ppanr\Documents\GitHub\Braille_kids')
addpath('C:\Users\ppanr\Documents\GitHub\Braille_kids\Violinplot-Matlab-master')
ft_defaults;

%% set directories and sub numbers
project_dir =  'R:\DRS-KidsOPM\Temp_warp_paper\pseudoMRI\';
project_dir_indiv = 'R:\DRS-KidsOPM\Temp_warp_paper\pseudomri_no_coreg_error\';
datadir = [project_dir,'Data',filesep,'BIDS',filesep];
datadir_indiv = [project_dir_indiv,'Data',filesep,'BIDS',filesep];

extension_T = '';
good_subs = [1:13 15 16 18 21 24:26];
f_GB = figure;f_GB.Color = 'w';set(f_GB,"Units","normalized","Position",[0,0,1,1])

path_mni_mesh = 'R:\DRS-KidsOPM\';
files_mni_mesh = 'meshes_mni.mat';
load([path_mni_mesh files_mni_mesh]); 

%% Loop through subjects, plotting peaks on glass brain
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.1,'edgecolor','none')
hold on
axis equal
view([-120 30])
fig = gcf;
fig.Color = 'w';
ax = gca;
ax.FontSize = 14;
for sub_i = good_subs
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'

    path_Tstat = [datadir,'derivatives',filesep,'Tstats',extension_T,filesep,'sub-',sub,filesep];
    path_Tstat_indiv = [datadir_indiv,'derivatives',filesep,'Tstats',extension_T,filesep,'sub-',sub,filesep];

    D1_tstat = ft_read_mri([path_Tstat 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_index.nii']);
    D1_tstat_indiv = ft_read_mri([path_Tstat_indiv 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_index.nii']);
    D4_tstat = ft_read_mri([path_Tstat 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_pinky.nii']);
    D4_tstat_indiv = ft_read_mri([path_Tstat_indiv 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_pinky.nii']);
    D1_peak(sub_i) = min(D1_tstat.anatomy(:));
    D1_peak_indiv(sub_i) = min(D1_tstat_indiv.anatomy(:));
    D4_peak(sub_i) = min(D4_tstat.anatomy(:));
    D4_peak_indiv(sub_i) = min(D4_tstat_indiv.anatomy(:));

    fid = fileread([path_Tstat,'D1_peak_mni.txt']);
    dip_loc_D1_mm = str2num(fid(43:end));
    fid = fileread([path_Tstat,'D4_peak_mni.txt']);
    dip_loc_D4_mm = str2num(fid(43:end));
    dip_loc_D1 = dip_loc_D1_mm./1000;
    dip_loc_D4 = dip_loc_D4_mm./1000;

    fid = fileread([path_Tstat_indiv,'D1_peak_mni.txt']);
    tmp = str2num(fid(43:end));
    dip_loc_D1_mm_indiv = tmp(1,:);
    fid = fileread([path_Tstat_indiv,'D4_peak_mni.txt']);
    tmp = str2num(fid(43:end));
    dip_loc_D4_mm_indiv = tmp(1,:);
    dip_loc_D1_indiv = dip_loc_D1_mm_indiv./1000;
    dip_loc_D4_indiv = dip_loc_D4_mm_indiv./1000;

    D1(sub_i,:) = dip_loc_D1;
    D1_indiv(sub_i,:) = dip_loc_D1_indiv;
    D4(sub_i,:) = dip_loc_D4;
    D4_indiv(sub_i,:) = dip_loc_D4_indiv;

    distance_D1(sub_i) = sqrt(sum((dip_loc_D1_mm-dip_loc_D1_mm_indiv).^2));
    distance_D4(sub_i) = sqrt(sum((dip_loc_D4_mm-dip_loc_D4_mm_indiv).^2));
    h_index = plot3(dip_loc_D1(1),dip_loc_D1(2),dip_loc_D1(3),'rx',...
        'markersize',10,'LineWidth',3,'DisplayName','D1');
    h_pinky = plot3(dip_loc_D4(1),dip_loc_D4(2),dip_loc_D4(3),'r+',...
        'markersize',10,'LineWidth',3,'DisplayName','D4');
    h_index_indiv = plot3(dip_loc_D1_indiv(1),dip_loc_D1_indiv(2),dip_loc_D1_indiv(3),'bx',...
        'markersize',10,'LineWidth',3,'DisplayName','D1');
    h_pinky_indiv = plot3(dip_loc_D4_indiv(1),dip_loc_D4_indiv(2),dip_loc_D4_indiv(3),'b+',...
        'markersize',10,'LineWidth',3,'DisplayName','D4');
%     legend
%     [k2,~] = convhull(sourcepos,'Simplify',true);
%     trisurf(k2,sourcepos(:,1),sourcepos(:,2),sourcepos(:,3),...
%         'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.5,'EdgeColor','none',...
%         'DisplayName','Sensorimotor VOI')
end

% Remove empty subjects
distance_D1(distance_D1==0)=[];
distance_D4(distance_D4==0)=[];
D1_peak_indiv(D1_peak_indiv==0)=[];
D1_peak(D1_peak==0)=[];
D4_peak_indiv(D4_peak_indiv==0)=[];
D4_peak(D4_peak==0)=[];
% Take mean
mean_distance_D1 = mean(distance_D1)
mean_distance_D4 = mean(distance_D4)
mean_D1 = mean(D1,1);
mean_D1_indiv = mean(D1_indiv,1);
mean_D4 = mean(D4,1);
mean_D4_indiv = mean(D4_indiv,1);

dist_D1 = sqrt(sum((mean_D1-mean_D1_indiv).^2))
dist_D4 = sqrt(sum((mean_D4-mean_D4_indiv).^2))

% Plot averages
figure;
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.1,'edgecolor','none')
hold on
axis equal
view([-120 30])
fig = gcf;
fig.Color = 'w';
ax = gca;
ax.FontSize = 14;
h_index = plot3(mean_D1(1),mean_D1(2),mean_D1(3),'rx',...
    'markersize',10,'LineWidth',3,'DisplayName','D1');
h_pinky = plot3(mean_D4(1),mean_D4(2),mean_D4(3),'r+',...
    'markersize',10,'LineWidth',3,'DisplayName','D4');
h_index_indiv = plot3(mean_D1_indiv(1),mean_D1_indiv(2),mean_D1_indiv(3),'bx',...
    'markersize',10,'LineWidth',3,'DisplayName','D1');
h_pinky_indiv = plot3(mean_D4_indiv(1),mean_D4_indiv(2),mean_D4_indiv(3),'b+',...
    'markersize',10,'LineWidth',3,'DisplayName','D4');

%% Plot vectors on glass brain from true anatomy to warped template anatomy
figure(174)
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.1,'edgecolor','none')
hold on
axis equal
view([-120 30])
fig = gcf;
fig.Color = 'w';
ax = gca;
ax.FontSize = 14;
vect_D1 = D1 - D1_indiv;
vect_D4 = D4 - D4_indiv;

quiver3(D1(:,1),D1(:,2),D1(:,3),vect_D1(:,1),vect_D1(:,2),vect_D1(:,3),'r','LineWidth',2);
quiver3(D4(:,1),D4(:,2),D4(:,3),vect_D4(:,1),vect_D4(:,2),vect_D4(:,3),'b','LineWidth',2);

%% Bland-Altman plots

distance_D1_demean = distance_D1 - mean(distance_D1);
distance_D4_demean = distance_D4 - mean(distance_D4);

peak_D1_sd = std(D1_peak_indiv-D1_peak);
peak_D4_sd = std(D4_peak_indiv-D4_peak);

figure;
subplot(1,2,1);
plot((D1_peak_indiv+D1_peak)./2,D1_peak_indiv-D1_peak,'ko');hold on;plot((D1_peak_indiv+D1_peak)./2,zeros(size([1:20])),'k-')
plot((D1_peak_indiv+D1_peak)./2,zeros(size([1:20]))+1.96*peak_D1_sd,'r-');plot((D1_peak_indiv+D1_peak)./2,zeros(size([1:20]))-1.96*peak_D1_sd,'r-')
title('Index');xlabel('Mean pseudo-T statistic');ylabel('Difference between pseudo-T statistics')

subplot(1,2,2);
plot((D4_peak_indiv+D4_peak)./2,D4_peak_indiv-D4_peak,'ko');hold on;plot((D4_peak_indiv+D4_peak)./2,zeros(size([1:20])),'k-')
plot((D4_peak_indiv+D4_peak)./2,zeros(size([1:20]))+1.96*peak_D4_sd,'r-');plot((D4_peak_indiv+D4_peak)./2,zeros(size([1:20]))-1.96*peak_D4_sd,'r-')
title('Little');xlabel('Mean pseudo-T statistic');ylabel('Difference between pseudo-T statistics')

%% Tstat v distance

figure;
subplot(1,2,1)
plot(D1_peak_indiv, distance_D1,'ko');
title('Index');xlabel('"Real" pseudo-T stat');ylabel('Distance between peaks')
ylim([0 70])
subplot(1,2,2)
plot(D4_peak_indiv, distance_D4,'ko');
title('Little');xlabel('"Real" pseudo-T stat');ylabel('Distance between peaks')
ylim([0 70])

%% Bland-Altman plot real vs mean

distance_D1_demean = distance_D1 - mean(distance_D1);
distance_D4_demean = distance_D4 - mean(distance_D4);

peak_D1_sd = std(D1_peak_indiv-D1_peak);
peak_D4_sd = std(D4_peak_indiv-D4_peak);
figure;
subplot(1,2,1);
plot(D1_peak_indiv,D1_peak_indiv-D1_peak,'ko');hold on;yline(mean(D1_peak_indiv-D1_peak),'k-')
yline(mean(D1_peak_indiv-D1_peak)+1.96*peak_D1_sd,'r-');yline(mean(D1_peak_indiv-D1_peak)-1.96*peak_D1_sd,'r-')
yline(0,'--');%xlim([min(D1_peak_indiv) max(D1_peak_indiv)])
title('Index');xlabel('"Real" pseudo-T statistic');ylabel('Difference between pseudo-T statistics')

subplot(1,2,2);
plot(D4_peak_indiv,D4_peak_indiv-D4_peak,'ko');hold on;yline(mean(D4_peak_indiv-D4_peak),'k-')
yline(mean(D4_peak_indiv-D4_peak)+1.96*peak_D4_sd,'r-');yline(mean(D4_peak_indiv-D4_peak)-1.96*peak_D4_sd,'r-')
yline(0,'--');%xlim([min(D4_peak_indiv) max(D4_peak_indiv)])
title('Little');xlabel('"Real" pseudo-T statistic');ylabel('Difference between pseudo-T statistics')

%% Violin plot
addpath('C:\Users\ppynr2\OneDrive - The University of Nottingham\matlab\tools\Violinplot-Matlab-master')

figure
violins = violinplot(distance_D1);
ylim([0 22]); ylabel('Distance between individual and pseudo-MRI') ;xlim([0.6 1.4])

figure;
boxplot(distance_D1); hold on; plot(ones(1,length(distance_D1)),distance_D1,'k.')
