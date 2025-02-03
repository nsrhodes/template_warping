% Plot connectivity results for template warping comparison
%% Housekeeping
clear all
clc
close all

addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
addpath('C:\Users\ppanr\Documents\GitHub\Braille_kids')
ft_defaults;

%% set directories and sub numbers
project_dir =  'R:\DRS-KidsOPM\Temp_warp_paper\pseudoMRI\';
project_dir_indiv = 'R:\DRS-KidsOPM\Temp_warp_paper\individual\';
datadir = [project_dir,'Data',filesep,'BIDS',filesep];
datadir_indiv = [project_dir_indiv,'Data',filesep,'BIDS',filesep];

good_subs = [1:13 15 16 18 21 24:26];

path_mni_mesh = 'R:\DRS-KidsOPM\';
files_mni_mesh = 'meshes_mni.mat';
load([path_mni_mesh files_mni_mesh]); 

%% Loop through subjects, plotting connectivity on glass brain
figure(111)
figure(222)
i = 0;
for sub_i = good_subs
    i = i+1; %iterate through subplots
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
    path_AEC = [datadir,'derivatives',filesep,'AEC',filesep,'sub-',sub,filesep];
    path_AEC_indiv = [datadir_indiv,'derivatives',filesep,'AEC',filesep,'sub-',sub,filesep];
    
    files_AEC = ['sub-' sub '_' 'ses-001_task-braille_run-001_AEC_13_30_Hz_Z.mat'];
    
    AEC = load([path_AEC files_AEC]);
    AEC_indiv = load([path_AEC_indiv files_AEC]);
    AEC_all(:,:,sub_i) = AEC.AEC;
    AEC_indiv_all(:,:,sub_i) = AEC_indiv.AEC;
    figure(111)
    subplot(4,5,i)
    go_netviewer_perctl(AEC.AEC,0.95) %top 5% of connections
    figure(222)
    subplot(4,5,i)
    go_netviewer_perctl(AEC_indiv.AEC,0.95)


end

% Average connectivity
AEC_mean = mean(AEC_all,3);
AEC_indiv_mean = mean(AEC_indiv_all,3);

% Plot averages
figure
subplot(1,2,1)
go_netviewer_perctl(AEC_mean,0.95)
subplot(1,2,2)
go_netviewer_perctl(AEC_indiv_mean,0.95)
figure
subplot(1,2,1)
go_netviewer_perctl(AEC_mean,0.95)
view([1 0 0])
subplot(1,2,2)
go_netviewer_perctl(AEC_indiv_mean,0.95)
view([1 0 0])
figure
subplot(1,2,1)
go_netviewer_perctl(AEC_mean,0.95)
view([0 1 0])
subplot(1,2,2)
go_netviewer_perctl(AEC_indiv_mean,0.95)
view([0 1 0])

figure('color','w'); subplot(1,2,1); imagesc(AEC_mean); axis square
yticks([5, 14, 25, 37, 44, 53, 64, 76]); 
yticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
    'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
xticks([5, 14, 25, 37, 44, 53, 64, 76]); 
xticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
    'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
xtickangle(45)
h = colorbar; ylabel(h,'AEC')
set(gca,'Fontsize',12,'LineWidth',2.5,'TickLength',[0.025,0.025]);
subplot(1,2,2); imagesc(AEC_indiv_mean); axis square
yticks([5, 14, 25, 37, 44, 53, 64, 76]); 
yticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
    'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
xticks([5, 14, 25, 37, 44, 53, 64, 76]); 
xticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
    'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
xtickangle(45)
h = colorbar; ylabel(h,'AEC')
set(gca,'Fontsize',12,'LineWidth',2.5,'TickLength',[0.025,0.025]);
% Plot relationship between global connectivity with template and
% individual anatomical
figure
plot(sum(reshape(AEC_indiv_all,[6084,26]),1,'omitnan')/2,sum(reshape(AEC_all,[6084,26]),1,'omitnan')/2,'ko')
xlabel('Individual MRI GC');ylabel('Pseudo-MRI GC')
axis([0 750 0 750])
hold on;plot(xlim,ylim,'k--')

%% 
% Vectorise the connectivity matrices
figure
nn=0;
for sub_i = good_subs
    nn = nn+1;
    connect_vect = [];
    connect_vect_indiv = [];
for i = 1:78
    connect_vect = [connect_vect squeeze(AEC_all(i,i+1:end,sub_i))];
    connect_vect_indiv = [connect_vect_indiv squeeze(AEC_indiv_all(i,i+1:end,sub_i))];
end
connect_vect_all(sub_i,:) = connect_vect;
connect_vect_all_indiv(sub_i,:) = connect_vect_indiv;
subplot(4,5,nn)
plot(connect_vect_indiv,connect_vect,'ko')
xlabel('Individual MRI');ylabel('Pseudo-MRI')
xlim([-0.1 0.4]);ylim([-0.1 0.4])
hold on;plot([-0.1:0.1:0.4],[-0.1:0.1:0.4],'k--')
end
% Correlate connectivity vectors
connect_vect_all([14 17 19 20 22 23],:) = [];
connect_vect_all_indiv([14 17 19 20 22 23],:) = [];

for i = 1:20
    [r,p] = corrcoef(connect_vect_all(i,:), connect_vect_all_indiv(i,:));
    corr_connect_all(i) = r(1,2);
end

% Plot relationship between global connectivity with template and
% individual anatomical
figure
plot(sum(connect_vect_all_indiv,2,'omitnan'),sum(connect_vect_all,2,'omitnan'),'ko')
xlabel('Individual MRI GC');ylabel('Pseudo-MRI GC')
axis([0 600 0 600])
hold on;plot(xlim,ylim,'k--')
axis square

AEC_mean_vect = triu(AEC_mean,1);
AEC_mean_vect(AEC_mean_vect==0)=[];
AEC_indiv_mean_vect = triu(AEC_indiv_mean,1);
AEC_indiv_mean_vect(AEC_indiv_mean_vect==0)=[];

[ra,pa] = corrcoef(AEC_indiv_mean_vect,AEC_mean_vect);

[r_gc,p_gc] = corrcoef(sum(connect_vect_all_indiv,2,'omitnan'),sum(connect_vect_all,2,'omitnan'))

figure
plot(connect_vect_all_indiv(:),connect_vect_all(:),'ko')
xlim([-0.1 0.7]);ylim([-0.1 0.7])
xlabel('Individual MRI connectivity');ylabel('Pseudo-MRI connectivity')
hold on;plot(xlim,ylim,'k--')

[r_allconn, p_allcon] = corrcoef(connect_vect_all_indiv(:),connect_vect_all(:))