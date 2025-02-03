% Convert 4mm brain and 4mm tstat units from m to mm for moving to MNI
% space
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

extension_T = '_4mm';
good_subs = [1:13 15 16 18 21 24:26];
path_mni_mesh = 'R:\DRS-KidsOPM\';
files_mni_mesh = 'meshes_mni.mat';
load([path_mni_mesh files_mni_mesh]); 

%% Loop through subjects loading in brain and tstats
count = 0;
for sub_i = good_subs
    count = count+1;
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
    % set paths
    path_Tstat = [datadir,'derivatives',filesep,'Tstats',extension_T,filesep,'sub-',sub,filesep];
    path_Tstat_indiv = [datadir_indiv,'derivatives',filesep,'Tstats',extension_T,filesep,'sub-',sub,filesep];

    path_brain = [datadir,'derivatives',filesep,'sourcespace',filesep,'sub-',sub,filesep];
    path_brain_indiv = [datadir_indiv,'derivatives',filesep,'sourcespace',filesep,'sub-',sub,filesep];
    
    % 4mm warped brain
    brain = ft_read_mri([path_brain 'sub-' sub '_brain_4mm.nii']);
    brain.transform(1:3,1:4) = brain.transform(1:3,1:4).*1e3;
    ft_write_mri([path_brain 'sub-' sub '_brain_4mm_mm.nii'],brain,'dataformat','nifti');
    % 4mm individual brain
    brain_indiv = ft_read_mri([path_brain_indiv 'sub-' sub '_brain_4mm.nii']);
    brain_indiv.transform(1:3,1:4) = brain_indiv.transform(1:3,1:4).*1e3;
    ft_write_mri([path_brain_indiv 'sub-' sub '_brain_4mm_mm.nii'],brain_indiv,'dataformat','nifti');
    % 4mm D1 warped tstat
    D1_tstat = ft_read_mri([path_Tstat 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_index.nii']);
    D1_tstat.transform(1:3,1:4) = D1_tstat.transform(1:3,1:4).*1e3;
    ft_write_mri([path_Tstat 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_index_mm.nii'],D1_tstat, 'dataformat','nifti');
    % 4mm D1 individual tstat
    D1_tstat_indiv = ft_read_mri([path_Tstat_indiv 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_index.nii']);
    D1_tstat_indiv.transform(1:3,1:4) = D1_tstat_indiv.transform(1:3,1:4).*1e3;
    ft_write_mri([path_Tstat_indiv 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_index_mm.nii'],D1_tstat_indiv,'dataformat','nifti');
    % 4mm D4 warped tstat
    D4_tstat = ft_read_mri([path_Tstat 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_pinky.nii']);
    D4_tstat.transform(1:3,1:4) = D4_tstat.transform(1:3,1:4).*1e3;
    ft_write_mri([path_Tstat 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_pinky_mm.nii'],D4_tstat,'dataformat','nifti');
    % 4mm D4 individal tstat
    D4_tstat_indiv = ft_read_mri([path_Tstat_indiv 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_pinky.nii']);
    D4_tstat_indiv.transform(1:3,1:4) = D4_tstat_indiv.transform(1:3,1:4).*1e3;
    ft_write_mri([path_Tstat_indiv 'sub-' sub '_ses-001_task-braille_run-001_pseudoT_pinky_mm.nii'],D4_tstat_indiv,'dataformat','nifti');
  
end

