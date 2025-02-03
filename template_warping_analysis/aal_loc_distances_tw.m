% Plot connectivity results for template warping comparison
%% Housekeeping
clear all
clc
close all

addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
ft_defaults;
addpath('C:\Users\ppanr\Documents\GitHub\Braille_kids\BrainPlots')
addpath('C:\Users\ppanr\Documents\GitHub\Braille_kids\gifti-1.8\')
%% set directories and sub numbers
project_dir =  'R:\DRS-KidsOPM\Temp_warp_paper\pseudoMRI\';
project_dir_indiv = 'R:\DRS-KidsOPM\Temp_warp_paper\pseudomri_no_coreg_error\';
datadir = [project_dir,'Data',filesep,'BIDS',filesep];
datadir_indiv = [project_dir_indiv,'Data',filesep,'BIDS',filesep];

good_subs = [1:13 15 16 18 21 24:26];

% path_mni_mesh = 'R:\DRS-KidsOPM\';
% files_mni_mesh = 'meshes_mni.mat';
% load([path_mni_mesh files_mni_mesh]); 

%% Get AAL centroids
i = 0;
figure
for sub_i = good_subs
    i = i+1 %iterate through subplots
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
    path_AEC = [datadir,'derivatives',filesep,'AEC',filesep,'sub-',sub,filesep];
    path_AEC_indiv = [datadir_indiv,'derivatives',filesep,'AEC',filesep,'sub-',sub,filesep];
    
    files_AEC = ['sub-' sub '_' 'ses-001_task-braille_run-001_lead_fields.mat'];
    
    load([path_AEC files_AEC]);
    all_sourcepos(:,:,i) = sourcepos;
    senspos{i} = S.sensor_info.pos;
    all_lf{i} = Lead_fields(:);
    load([path_AEC_indiv files_AEC]);
    all_sourcepos_indiv(:,:,i) = sourcepos;
    senspos_indiv{i} = S.sensor_info.pos;
    all_lf_indiv{i} = Lead_fields(:);

    %get rot sensors
%     [R(:,:,i),t(:,i)] = get_rot(senspos_indiv{i},senspos{i});
%     senspos_indiv_trans = (R(:,:,i)*senspos_indiv{i}'+t(:,i))';
%     all_sourcepos_indiv(:,:,i) = (R(:,:,i)*all_sourcepos_indiv(:,:,i)'+t(:,i))';
%     hold on;
%     plot3(squeeze(all_sourcepos_indiv(:,1,i)),squeeze(all_sourcepos_indiv(:,2,i)),squeeze(all_sourcepos_indiv(:,3,i)),'o')
%     plot3(squeeze(all_sourcepos(:,1,i)),squeeze(all_sourcepos(:,2,i)),squeeze(all_sourcepos(:,3,i)),'x')
% 
%     plot3(senspos{i}(:,1),senspos{i}(:,2),senspos{i}(:,3),'kx','MarkerSize',10)
%     plot3(senspos_indiv_trans(:,1),senspos_indiv_trans(:,2),senspos_indiv_trans(:,3),'rx','MarkerSize',10)
% 
%     title(sub);
%     drawnow
%     pause(0.5)
end

for i = 1:20
    for j = 1:78
        sourcepos_dist(j,i) = sqrt(sum((all_sourcepos_indiv(j,:,i)-all_sourcepos(j,:,i)).^2));
    end
    [r_lf(i),p_lf(i)] = corr(all_lf{i},all_lf_indiv{i});
end

mean_sourcepos_dist = mean(sourcepos_dist,1);
mean_sourcepos_dist_mm = mean_sourcepos_dist*1e3;
all_mean = mean(mean_sourcepos_dist_mm);
all_se = std(mean_sourcepos_dist_mm)./sqrt(20);

all_sourcepos_mean = mean(all_sourcepos,3);
all_sourcepos_indiv_mean = mean(all_sourcepos_indiv,3);

all_sourcepos_mean_dist = mean(sqrt(sum((all_sourcepos_indiv_mean'-all_sourcepos_mean').^2)))*1000;


%% plot best and worst participant locations

% load in template warped mesh
[~,i_max] = max(mean_sourcepos_dist_mm);
[v,i_min] = min(mean_sourcepos_dist_mm);
i_max=3;
sub_i = good_subs(i_max);
sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
path_meshes = [datadir,'derivatives',filesep,'sourcespace',filesep,'sub-',sub,filesep];
files_meshes = ['sub-',sub,'_meshes.mat'];
meshes_max = load([path_meshes files_meshes],'meshes');

sub_i = good_subs(i_min);
sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
path_meshes = [datadir,'derivatives',filesep,'sourcespace',filesep,'sub-',sub,filesep];
files_meshes = ['sub-',sub,'_meshes.mat'];
meshes_min = load([path_meshes files_meshes],'meshes');

figure;
% max distance
% plot warped locations
subplot(1,2,1)
ft_plot_mesh(meshes_max.meshes,'facecolor',[.5 .5 .5],'facealpha',.1,'edgecolor','none');
hold on
scatter3(squeeze(all_sourcepos(:,1,i_max)),squeeze(all_sourcepos(:,2,i_max)),squeeze(all_sourcepos(:,3,i_max)),'red','filled')

% plot individual locations
scatter3(squeeze(all_sourcepos_indiv(:,1,i_max)),squeeze(all_sourcepos_indiv(:,2,i_max)),squeeze(all_sourcepos_indiv(:,3,i_max)),'blue','filled');

% min distance
% plot warped locations
subplot(1,2,2)
ft_plot_mesh(meshes_min.meshes,'facecolor',[.5 .5 .5],'facealpha',.1,'edgecolor','none');
hold on
plot3(squeeze(all_sourcepos(:,1,i_min)),squeeze(all_sourcepos(:,2,i_min)),squeeze(all_sourcepos(:,3,i_min)),'r.','MarkerSize',15)

% plot individual locations
plot3(squeeze(all_sourcepos_indiv(:,1,i_min)),squeeze(all_sourcepos_indiv(:,2,i_min)),squeeze(all_sourcepos_indiv(:,3,i_min)),'b.','MarkerSize',15);

%% Euclidean distance by region

addpath('C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\Braille_kids\BrainPlots')
addpath('C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\Braille_kids\gifti-1.8\')
figure; subplot(1,3,1)
PaintBrodmannAreas_chooseview(mean(sourcepos_dist*1e3,2), 78, 78, [1 5], 'Pseudo-T',[],[90 0 0])

subplot(1,3,2)
PaintBrodmannAreas_chooseview(mean(sourcepos_dist*1e3,2), 78, 78, [1 5], 'Pseudo-T',[],[0 90 0])
subplot(1,3,3)
PaintBrodmannAreas_chooseview(mean(sourcepos_dist*1e3,2), 78, 78, [1 5], 'Pseudo-T',[],[0 0 90])

mean_byreg = mean(mean(sourcepos_dist,2))*1e3
sd_byreg = std(mean(sourcepos_dist,2))*1e3

