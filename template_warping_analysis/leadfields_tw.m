% Plot connectivity results for template warping comparison
%% Housekeeping
clear all
clc
close all

addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
ft_defaults;

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
%% Save out AAL volumes

% i = 0;
% for sub_i = good_subs
%     i = i+1;
%     sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
%     path_meshes = [datadir,'derivatives',filesep,'sourcespace',filesep,'sub-',sub,filesep];
%     path_meshes_indiv = [datadir_indiv,'derivatives',filesep,'sourcespace',filesep,'sub-',sub,filesep];
% 
%     files_AAL_regions = ['sub-',sub,'_AAL_regions.nii.gz'];
%     AAL_regions = ft_read_mri([path_meshes,files_AAL_regions]);
%     for j = 1:78
%         AAL_vols(i,j) = length(find(AAL_regions.anatomy(:,:,:,j)>0));
%     end
%     AAL_vols_sub = AAL_vols(i,:);
%     save([path_meshes 'sub-' sub '_AAL_vols.mat'],'AAL_vols_sub')
%     clear AAL_regions AAL_vols_sub
%     AAL_regions = ft_read_mri([path_meshes_indiv,files_AAL_regions]);
%     for j = 1:78
%         AAL_vols_indiv(i,j) = length(find(AAL_regions.anatomy(:,:,:,j)>0));
%     end 
%     AAL_vols_sub = AAL_vols_indiv(i,:);
%     save([path_meshes_indiv 'sub-' sub '_AAL_vols.mat'],'AAL_vols_sub')
%     clear AAL_regions AAL_vols_sub
%     save([path_meshes 'sub-' sub '_AAL_vols.mat'],'AAL_vols')
% end

%% Read in AAL volume sizes

i = 0;
figure;
for sub_i = good_subs
    i = i+1;
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
    path_meshes = [datadir,'derivatives',filesep,'sourcespace',filesep,'sub-',sub,filesep];
    path_meshes_indiv = [datadir_indiv,'derivatives',filesep,'sourcespace',filesep,'sub-',sub,filesep];

    load([path_meshes 'sub-' sub '_AAL_vols.mat'])
    AAL_vols_tw(i,:) = AAL_vols(i,:);
    
    load([path_meshes_indiv 'sub-' sub '_AAL_vols.mat'])
    AAL_vols_indiv(i,:) = AAL_vols_sub;
    plot(AAL_vols_indiv(i,:),AAL_vols_tw(i,:),'ko');hold on;
end
[r,p] = corr(AAL_vols_indiv(:),AAL_vols_tw(:))

%% plot best and worst participant locations

% load in template warped mesh
[~,i_max] = max(mean_sourcepos_dist_mm);
[v,i_min] = min(mean_sourcepos_dist_mm);
i_max=1;
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
plot3(squeeze(all_sourcepos(:,1,i_max)),squeeze(all_sourcepos(:,2,i_max)),squeeze(all_sourcepos(:,3,i_max)),'r.','MarkerSize',15)

% plot individual locations
plot3(squeeze(all_sourcepos_indiv(:,1,i_max)),squeeze(all_sourcepos_indiv(:,2,i_max)),squeeze(all_sourcepos_indiv(:,3,i_max)),'b.','MarkerSize',15);

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
PaintBrodmannAreas_chooseview(mean(sourcepos_dist,2), 78, 78, [0 0.02], 'Pseudo-T',[],[90 0 0])

subplot(1,3,2)
PaintBrodmannAreas_chooseview(mean(sourcepos_dist,2), 78, 78, [0 0.02], 'Pseudo-T',[],[0 90 0])
subplot(1,3,3)
PaintBrodmannAreas_chooseview(mean(sourcepos_dist,2), 78, 78, [0 0.02], 'Pseudo-T',[],[0 0 90])

mean_byreg = mean(mean(sourcepos_dist,2))*1e3
sd_byreg = std(mean(sourcepos_dist,2))*1e3
%% Load in peak distances
extension_T = '';
i = 0;
for sub_i = good_subs
    i = i+1;
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'

    path_Tstat = [datadir,'derivatives',filesep,'Tstats',extension_T,filesep,'sub-',sub,filesep];
    path_Tstat_indiv = [datadir_indiv,'derivatives',filesep,'Tstats',extension_T,filesep,'sub-',sub,filesep];

    fid = fileread([path_Tstat,'D1_peak_mni.txt']);
    dip_loc_D1_mm = str2num(fid(43:end));
    dip_loc_D1 = dip_loc_D1_mm./1000;

    fid = fileread([path_Tstat_indiv,'D1_peak_mni.txt']);
    tmp = str2num(fid(43:end));
    dip_loc_D1_mm_indiv = tmp(1,:);
    dip_loc_D1_indiv = dip_loc_D1_mm_indiv./1000;

    D1(i,:) = dip_loc_D1;
    D1_indiv(i,:) = dip_loc_D1_indiv;
    distance_D1(i) = sqrt(sum((dip_loc_D1_mm-dip_loc_D1_mm_indiv).^2));
end
%% Load in connectivity
i=0;
for sub_i = good_subs
    i = i+1 %iterate through subplots
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
    path_AEC = [datadir,'derivatives',filesep,'AEC',filesep,'sub-',sub,filesep];
    path_AEC_indiv = [datadir_indiv,'derivatives',filesep,'AEC',filesep,'sub-',sub,filesep];
    
    files_AEC = ['sub-' sub '_' 'ses-001_task-braille_run-001_AEC_13_30_Hz_Z.mat'];
    
    AEC = load([path_AEC files_AEC]);
    AEC_indiv = load([path_AEC_indiv files_AEC]);
    AEC_all(:,:,i) = AEC.AEC;
    AEC_indiv_all(:,:,i) = AEC_indiv.AEC;
end

% Average connectivity
AEC_mean = mean(AEC_all,3);
AEC_indiv_mean = mean(AEC_indiv_all,3);

%% 
% Vectorise the connectivity matrices
i=0;
for sub_i = good_subs
    i = i+1;
    this_vect = [];this_vect_indiv=[];
    connect_vect = AEC_all(:,:,i);
    connect_vect_indiv = AEC_indiv_all(:,:,i);
    this_vect = connect_vect(:);
    this_vect_indiv = connect_vect_indiv(:);
    nan_inds = find(isnan(this_vect));
    this_vect(nan_inds) = [];
    this_vect_indiv(nan_inds) = [];
    connect_vect_all(i,:) = this_vect;
    connect_vect_all_indiv(i,:) = this_vect_indiv;

end

% Correlate connectivity vectors
for i = 1:20
connectivity_corr(i) = corr(connect_vect_all(i,:)',connect_vect_all_indiv(i,:)');
end 
%%

f1 = figure;
fwidth = 350;
fheight = 300;
f1.Position([3,4]) = [fwidth,fheight];
plot(mean_sourcepos_dist_mm,distance_D1,'ko');
xlabel('Distance between AAL regions (mm)');
ylabel('Distance between peaks (mm)')
[r_struc2peak,p_struc2peak]= corr(mean_sourcepos_dist_mm',distance_D1')
set(gca, 'FontSize',12)
X = [ones(length(mean_sourcepos_dist_mm'),1) mean_sourcepos_dist_mm'];
y = distance_D1';
b = X\y;
refl = refline(b(2),b(1));
axis square
box off

f2 = figure;
f2.Position([3,4]) = [fwidth,fheight];
plot(r_lf,distance_D1,'ko');
xlabel('Lead field correlation');
ylabel('Distance between peaks (mm)')
[r_lf2peak,p_lf2peak]= corr(r_lf',distance_D1')
set(gca, 'FontSize',12)
X = [ones(length(r_lf'),1) r_lf'];
y = distance_D1';
b = X\y;
ref2 = refline(b(2),b(1));
axis square
box off

f3 = figure;
f3.Position([3,4]) = [fwidth,fheight];
plot(mean_sourcepos_dist_mm,connectivity_corr,'ko');
xlabel('Distance between AAL regions (mm)');
ylabel('Connectivity correlation')
set(gca, 'FontSize',12)
X = [ones(length(mean_sourcepos_dist_mm'),1) mean_sourcepos_dist_mm'];
y = connectivity_corr';
b = X\y;
ref3 = refline(b(2),b(1));
axis square
box off
[r_struc2conn,p_struc2conn] = corr(mean_sourcepos_dist_mm',connectivity_corr')

f4 = figure;
f4.Position([3,4]) = [fwidth,fheight];
plot(r_lf,connectivity_corr,'ko');
xlabel('Lead field correlation');
ylabel('Connectivity correlation')
set(gca, 'FontSize',12)
X = [ones(length(r_lf'),1) r_lf'];
b = X\y;
ref4 = refline(b(2),b(1));
axis square
box off
[r_lf2conn,p_lf2conn] = corr(r_lf',connectivity_corr')

%% brain volume by AAL region
figure; subplot(1,3,1)
PaintBrodmannAreas_chooseview(mean(AAL_vols_indiv,1), 78, 78, [0 50000], 'Pseudo-T',[],[90 0 0])

subplot(1,3,2)
PaintBrodmannAreas_chooseview(mean(AAL_vols_indiv,1), 78, 78, [0 50000], 'Pseudo-T',[],[0 90 0])
subplot(1,3,3)
PaintBrodmannAreas_chooseview(mean(AAL_vols_indiv,1), 78, 78, [0 50000], 'Pseudo-T',[],[0 0 90])


