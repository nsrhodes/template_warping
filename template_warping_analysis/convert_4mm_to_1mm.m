%% Housekeeping
clear all
clc
close all

addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
ft_defaults;


%%
path = 'R:\DRS-KidsOPM\Temp_warp_paper\individual\Data\BIDS\derivatives\Tstats_4mm\sub-103\sub-103_ses-001_task-braille_run-001_pseudoT_pinky.nii'
brain = ft_read_mri('R:\DRS-KidsOPM\Temp_warp_paper\individual\Data\BIDS\derivatives\sourcespace\sub-103\sub-103_brain.nii')
tstat = ft_read_mri(path);
tstat.transform(1:3,1:4) = tstat.transform(1:3,1:4).*1e3;
tstat = ft_convert_units(tstat,'m');
x = tstat.transform(1,4) + (0:size(tstat.anatomy,1)-1) * tstat.hdr.xsize;
y = tstat.transform(2,4) + (0:size(tstat.anatomy,2)-1) * tstat.hdr.ysize;
z = tstat.transform(3,4) + (0:size(tstat.anatomy,3)-1) * tstat.hdr.zsize;

cfg = [];
cfg.xrange = [min(x) max(x)];
cfg.yrange = [min(y) max(y)];
cfg.zrange = [min(z) max(z)];
%cfg.dim = [256 256 256];
cfg.resolution = 0.001;
cfg.method = 'linear';
[upsample_tstat] = ft_volumereslice(cfg, tstat);
upsample_tstat = ft_convert_units(upsample_tstat,'mm');
% upsample_tstat.transform = tstat.transform*1e3;
% upsample_tstat.transform(1,1) = upsample_tstat.transform(1,1)/4;
% upsample_tstat.transform(2,2) = upsample_tstat.transform(2,2)/4;
% upsample_tstat.transform(3,3) = upsample_tstat.transform(3,3)/4;
% upsample_tstat.transform(4,4) = 1;
ft_write_mri('R:\DRS-KidsOPM\Temp_warp_paper\individual\Data\BIDS\derivatives\Tstats_4mm\sub-103\Tstat_4mm_in_1mm_pinky.nii',upsample_tstat,...
    'dataformat','nifti');