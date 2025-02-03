%% Convert einscan of headshape into a filled-in image for flirt coregistration to template anatomical MRI %%
%% Housekeeping

close all
clear all
clc


project_dir = 'R:\DRS-KidsOPM\Paediatric_OPM_Notts_AdultData\';
wsl_project_dir = '/mnt/r/DRS-KidsOPM/Paediatric_OPM_Notts_AdultData/';


% Setup FSL directory via wsl
wsl_fsl_dir = '/home/ppynr2/fsl';
if system(['wsl cd ',wsl_fsl_dir]) ~= 0
    error('Wrong fsl directory')
end

% Fieldtrip
addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
ft_defaults

% Get original template mesh - only needed if MRI has not been used before
template_mesh = 1;
convert_einscan_units = 1;

% pad template MRI for space to warp to larger heads without clipping
padded = 1;
% optional crop template MRI to remove chin etc. for better fitting
crop = 0;

%% Get original template mesh

if template_mesh
    cd([project_dir,'Data',filesep,'BIDS',filesep,'derivatives',filesep,'template_warping']) % cd to directory
    [mri_file,mri_path] = uigetfile('*.nii'); % pick MRI nifty
    mri_orig = ft_read_mri([mri_path mri_file]); % read in MRI
    if padded 
    crg = [];
    cfg.dim = mri_orig.dim + 50;
    cfg.method = 'cubic';
    cfg.resolution = 1;
    mri_orig = ft_volumereslice(cfg,mri);
    ft_write_mri([mri_path mri_file(1:end-4) '_padded.nii'],mri_new,'dataformat','nifti')
    end 
    cfg = [];
    cfg.output    = {'brain' 'scalp'};
    cfg.scalpsmooth = 10;
    %cfg.scalpthreshold = 0.03;
    segmentedmri  = ft_volumesegment(cfg, mri_orig); % segment MRI
    disp('Done!')
    cfg = [];cfg.tissue = {'brain' 'scalp'};
    cfg.numvertices = [5000 5000 5000];
    mesh2 = ft_prepare_mesh(cfg,segmentedmri); % get meshes
    mesh1 = ft_convert_units(mesh2,'m');
    for n = 1:size(mesh1,2)
        meshes(n).pnt = mesh1(n).pos;
        meshes(n).tri = mesh1(n).tri;
        meshes(n).unit = mesh1(n).unit;
        meshes(n).name = cfg.tissue{n};
    end

    % plot meshes
    figure;ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
    % turn into pointcloud
    plypts = meshes(2).pnt;
    plycolour = repmat([0.5 0.5 0.5],length(plypts),1);
    meshptc = pointCloud(plypts,'Color',plycolour);
    pcshow(meshptc)
    % Write out scalp points
    pcwrite(meshptc,[mri_path mri_file(1:8) 'template_points.ply'],'PLYformat','binary');
end

%% Convert einscan units
% if directlty from einscan, units need to converted to m
if convert_einscan_units
    cd([project_dir,'Data',filesep,'BIDS',filesep,'derivatives',filesep,'einscans']) % cd to directory
    [fname,path] = uigetfile('.ply'); % select einscan of head
    sf = 0.001; % scale factor is converting from mm to m

    einscan = pcread([path filesep fname]);
    einscan_new = pointCloud(einscan.Location*sf);
    einscan_new.Color = einscan.Color;
    % Write out einscan with converted units
    pcwrite(einscan_new,[path filesep fname(1:end-4) '_converted.ply'])
end
% GO TO MESHLAB AND TRANSFORM CONVERTED EINSCAN TO TEMPLATE SCALP POINTS
% AND SAVE OUT TRANSFORMATION FOR NEXT STEPS

%% Load in einscan

cd([project_dir,'Data',filesep,'BIDS',filesep,'derivatives',filesep,'template_warping']) % cd to directory
[filename,einscanpath] = uigetfile('*.ply'); % select cropped einscan image
sf = 1; % scalefactor to convert to metres
einscan = pcread([einscanpath filename]); % load in einscan
einscan_new = pointCloud(einscan.Location*sf); % convert to m
einscan_new.Color = einscan.Color;
h = figure(102);showPointCloud(einscan_new); % show einscan
einscan_points = double(einscan_new.Location); % extract points

%% Apply crude transform from einscan to template from meshlab

einscan2template = load([einscanpath filename(1:7) '_einscan2template.txt']); % read in transformation text file
einscan_points = [einscan2template(1:3,1:3)*einscan_points']' + einscan2template(1:3,4)'; % apply transform
einscan_points_mm = round(einscan_points.*1e3);
% convert einscan points to mm and centre einscan points
einscan_points_cent = mean(einscan_points_mm,1);

%% Load in template points

if ~exist('segmentedmri','var')
    cd([project_dir,'Data',filesep,'BIDS',filesep,'derivatives',filesep,'template_warping']) % cd to directory
    [mri_file,mri_path] = uigetfile('*.nii'); % pick MRI nifty
    mri_orig = ft_read_mri([mri_path mri_file]); % read in MRI
    cfg = [];
    cfg.output    = {'brain' 'scalp'};
    cfg.scalpsmooth = 10;
    %cfg.scalpthreshold = 0.03;
    segmentedmri  = ft_volumesegment(cfg, mri_orig); % segment MRI
end

%% Load in template image

if ~exist('mri_orig','var')
cd(einscanpath) % cd to directory
[mrifilename,~] = uigetfile('*.nii'); % select unwarped template file to get dimensions
mri = ft_read_mri(mrifilename);
else 
    mri = mri_orig;
end 

im_dim = mri.dim;
% centre einscan
einscan_points_mm_cent = einscan_points_mm - round(einscan_points_cent)+round(im_dim/2);

%% convert to binary images

einscan_image = zeros(im_dim);
einscan_points_mm_cent(einscan_points_mm_cent<=0) = 1;

for j = 1:size(einscan_points,1)
    if einscan_points_mm_cent(j,:) <= im_dim
        einscan_image(einscan_points_mm_cent(j,1),einscan_points_mm_cent(j,2),einscan_points_mm_cent(j,3)) = 1;
    end
end
% fill in using convex hull
einscan_image_filled = [];
for i = 1:im_dim(3)
    einscan_image_filled(:,:,i) = bwconvhull(squeeze(einscan_image(:,:,i)));
end
% define meshgrid for patch
[px,py,pz] = meshgrid(1:im_dim(2),1:im_dim(1),1:im_dim(3));
% plot patch to inspect
figure
pnew=patch(isosurface(px,py,pz,einscan_image_filled,0));
set(pnew, 'FaceColor', 'red', 'EdgeColor', 'none');
daspect([1 1 1])
view(3)
hold on
camlight
% save image
einscan_mri = mri;
einscan_mri.anatomy = einscan_image_filled;
ft_write_mri([einscanpath filename(1:7) '_einscan_img.nii'],einscan_mri,'dataformat','nifti');

% pad segmented mri rather than resegment
template_image = segmentedmri.scalp;
% fill in template scalp
template_image_filled = [];
for i = 1:im_dim(3)
    template_image_filled(:,:,i) = bwconvhull(squeeze(template_image(:,:,i)));
end
if crop
    % Define line to crop
    for i = 1:im_dim(2)
        template_image_filled(:,i,1:ceil(-0.33*i+92))=0;
    end
end
% plot patch to inspect
figure
pnew=patch(isosurface(px,py,pz,template_image_filled,0));
set(pnew, 'FaceColor', 'red', 'EdgeColor', 'none');
daspect([1 1 1])
view(3)
hold on
camlight
% save image
template_mri = mri;
template_mri.anatomy = double(template_image_filled);
ft_write_mri([einscanpath filename(1:7) '_template_img.nii'],template_mri,'dataformat','nifti');
%pause

%% FLIRT

sub = filename(1:7);
% write temporary file with fsl commands
f = fopen('tmp.sh','w');
fprintf(f,'#!/bin/bash\nexport FSLOUTPUTTYPE=NIFTI_GZ\n')
% first command - using FLIRT to warp template scalp to einscan scalp
str1 = [wsl_fsl_dir,'/bin/flirt -in ',wsl_project_dir,'Data/BIDS/derivatives/template_warping/' sub '/' sub '_template_img.nii -ref ',wsl_project_dir,'Data/BIDS/derivatives/template_warping/' sub '/' sub '_einscan_img.nii -out ',wsl_project_dir,'Data/BIDS/derivatives/template_warping/' sub '/' sub '_template_img_crg.nii -omat ',wsl_project_dir,'Data/BIDS/derivatives/template_warping/' sub '/' sub '_template_img_crg.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 -interp trilinear'];
% second command - applying transformation from first command to the padded
% template mri
if padded
    str2 = [wsl_fsl_dir,'/bin/flirt -in ',wsl_project_dir,'Data/BIDS/derivatives/template_warping/' sub '/' sub '_template_file_padded.nii -applyxfm -init ',wsl_project_dir,'Data/BIDS/derivatives/template_warping/' sub '/' sub '_template_img_crg.mat -out ',wsl_project_dir,'Data/BIDS/derivatives/template_warping/' sub '/' sub '_template_file_crg.nii -paddingsize 0.0 -interp trilinear -ref ',wsl_project_dir,'Data/BIDS/derivatives/template_warping/' sub '/' sub '_template_file_padded.nii'];
else
    str2 = [wsl_fsl_dir,'/bin/flirt -in ',wsl_project_dir,'Data/BIDS/derivatives/template_warping/' sub '/' sub '_template_file.nii -applyxfm -init ',wsl_project_dir,'Data/BIDS/derivatives/template_warping/' sub '/' sub '_template_img_crg.mat -out ',wsl_project_dir,'Data/BIDS/derivatives/template_warping/' sub '/' sub '_template_file_crg.nii -paddingsize 0.0 -interp trilinear -ref ',wsl_project_dir,'Data/BIDS/derivatives/template_warping/' sub '/' sub '_template_file.nii'];
end
fprintf(f,'%s\n',str1); %add for str2
fprintf(f,'%s\n',str2);
fclose(f) ;% if ans = 0, delete temp file
wsl_stdout = system('wsl ./tmp.sh');
if wsl_stdout == 0
    delete tmp.sh
else
    error('Error using flirt in WSL')
end
% Unzip warped mri
gunzip([project_dir,'Data',filesep,'BIDS',filesep,'derivatives',filesep,'template_warping', filesep, sub, filesep, sub, '_template_file_crg.nii.gz'])

%% AFTER FLIRT

% segment warped mri
cd(einscanpath)
[mriwarpedname,~] = uigetfile('*.nii');
mri_warped = ft_read_mri(mriwarpedname);
cfg = [];
cfg.output    = {'brain' 'scalp'};
cfg.scalpsmooth = 10;
segmentedmri  = ft_volumesegment(cfg, mri_warped);
disp('Done!')
cfg = [];cfg.tissue = {'brain' 'scalp'};
cfg.numvertices = [5000 5000 5000];
mesh2 = ft_prepare_mesh(cfg,segmentedmri);
mesh1 = ft_convert_units(mesh2,'m');
for n = 1:size(mesh1,2)
    meshes(n).pnt = mesh1(n).pos;
    meshes(n).tri = mesh1(n).tri;
    meshes(n).unit = mesh1(n).unit;
    meshes(n).name = cfg.tissue{n};
end
figure;ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
cd(einscanpath)
save('meshes.mat','meshes');save('segmentedmri.mat','segmentedmri');

plypts = meshes(2).pnt;
plycolour = repmat([0.5 0.5 0.5],length(plypts),1);
meshptc = pointCloud(plypts,'Color',plycolour);

% extract scalp points and save for coregistration
pcshow(meshptc)
pcwrite(meshptc,[filename(1:7) '_warped_scalp_points'],'PLYformat','binary');
