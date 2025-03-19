%% Convert einscan of headshape into a filled-in image for flirt coregistration to template anatomical MRI %%

close all
clear all
clc

%% Set paths and options here!
% Fieldtrip
fieldtrip_path = "path/to/fieldtrip-20250114";
% directory containing template MRI and Einscan .ply files
project_dir = './example/';

% pad template MRI for space to warp to larger heads without clipping
padded = 1;
% optional crop template MRI to remove chin etc. for better fitting
crop = 1;

%%
if ~exist(fieldtrip_path,"dir")
    error("Set FieldTrip path first!")
end

if ~exist(project_dir,"dir")
    error("project_dir not set correctly")
end

addpath(fieldtrip_path)
ft_defaults


%% Get original template mesh

cd(project_dir) % cd to directory
[mri_file,mri_path] = uigetfile({'*.nii;*.nii.gz'},"Pick the template MRI"); % pick MRI nifty
mri_file_path = [mri_path mri_file];
if endsWith(mri_file_path,'gz')
    extension_len = 7;
else
    extension_len = 4;
end
mri_f_head = mri_file(1:end-extension_len);
template_mesh_stl = [mri_path mri_f_head '_template_points.stl'];

mri_orig = ft_read_mri(mri_file_path); % read in MRI

if ~exist(template_mesh_stl,'file')

    if padded
        padded_mri = [mri_path mri_f_head '_padded.nii'];
        crg = [];
        cfg.dim = mri_orig.dim + 50;
        cfg.method = 'cubic';
        cfg.resolution = 1;
        mri_orig = ft_volumereslice(cfg,mri_orig);
        ft_write_mri(padded_mri,mri_orig,'dataformat','nifti')
        padded_str = 'padded';
    end

    if ~exist([mri_path mri_f_head padded_str '_segmentation.mat'],"file")
        cfg = [];
        cfg.output    = {'brain' 'scalp'};
        cfg.scalpsmooth = 10;
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
        save([mri_path mri_f_head padded_str '_segmentation.mat'],"meshes","segmentedmri")
    end
end
load([mri_path mri_f_head padded_str '_segmentation.mat'],"meshes","segmentedmri")

% plot meshes
figure;
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
% turn into pointcloud
v = meshes(strcmp({meshes.name},'scalp')).pnt;
f = meshes(strcmp({meshes.name},'scalp')).tri;
% Write out scalp points
C = ones(size(f)).*0.5;C(:,2:3) = 0;
p = patch('Faces',f,'Vertices',v,'FaceVertexCData',C,'FaceAlpha',0.1);
shading flat
% shading interp
stlwrite(template_mesh_stl,f,v)


%%
waitfor(msgbox(["* Open Meshlab";...
    "1. Import the template mesh (...)_template_points.stl";...
    " ";...
    "2. Import the head-only einscan";...
    "    If the tempate mesh disappears, you need to rescale the head-only mesh ";...
    "    by 0.001 (from mm to m)";...
    "    If any of the meshes appear dark, with little surface texture, ";...
    "    invert the faces' orientation by selecting";...
    "     Filters>Normals,Curvature and Orientation>Invert Faces Orientation";...
    " ";...
    "3. Roughly align the head-only mesh using 'Points-based gluing'";...
    "    Be careful to select the MRI mesh, then 'Glue here' to bring everything";...
    "    into MRI space.";...
    " ";...
    "4. Freeze the transform by selecting the head only mesh, right-clicking and ";...
    "  selecting 'Matrix: Freeze current matrix'";...
    " ";...
    "5. Crop the head-only mesh to remove points below the chin, to achieve a ";...
    "  similar shape to the template MRI Mesh. Removing points below the plane";...
    "  connecting the mouth and the ears can be advantageous."
    " ";...
    "6. Export the mesh with a sensible suffix (e.g. sub-101_head_cropped.ply)";...
    " ";...
    "Once done, click OK below!";],...
    "Rough Alignment and cropping in Meshlab"));
    
%% Load in einscan

[filename,einscanpath] = uigetfile('*.ply',"Choose cropped head-only Einscan."); % select cropped einscan image
einscan = pcread([einscanpath filename]); % load in einscan
h = figure("Name","Cropped Head Mesh");showPointCloud(einscan); % show einscan
einscan_points = double(einscan.Location); % extract points

einscan_points_mm = round(einscan_points.*1e3);
% convert einscan points to mm and centre einscan points
einscan_points_cent = mean(einscan_points_mm,1);

%% Load in template points

if ~exist('segmentedmri','var')
    if padded;pad_str = "padded";else;pad_str="";end
    [mri_file,mri_path] = uigetfile('*.nii',sprintf("Pick the %s template MRI",pad_str)); % pick MRI nifty
    mri_orig = ft_read_mri([mri_path mri_file]); % read in MRI
    cfg = [];
    cfg.output    = {'brain' 'scalp'};
    cfg.scalpsmooth = 10;
    segmentedmri  = ft_volumesegment(cfg, mri_orig); % segment MRI
end

%% Load in template image

if ~exist('mri_orig','var') || padded
    cd(einscanpath) % cd to directory
    [mrifilename,~] = uigetfile('*.nii',"Pick the unwarped template MRI (padded if using)"); % select unwarped template file to get dimensions
    mri = ft_read_mri(mrifilename);
else
    mri = mri_orig;
end

im_dim = mri.dim;
% bring einscan into voxel volume
einscan_points_mm_cent=round(ft_warp_apply(inv(mri.transform),einscan_points_mm));
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
einscan_filled_nii = [einscanpath mri_f_head '_einscan_img.nii'];
ft_write_mri(einscan_filled_nii,einscan_mri,'dataformat','nifti');

% pad segmented mri rather than resegment
template_image = segmentedmri.scalp;


%% Neck cropping
fig = uifigure("Name","Crop neck below the red line. Then click Confirm.");
ax1 = uiaxes(fig, 'Units','Normalized', 'Position',[.1, .3, .8, .7]);

overlay_image = squeeze(template_image(floor(size(template_image,1)/2),:,:));
overlay_image= overlay_image + squeeze(einscan_mri.anatomy(floor(size(einscan_mri.anatomy,1)/2),:,:));

imagesc(ax1,imrotate(overlay_image,-270))
hold(ax1,'on')
crop_intercept = 185.16;
if ~padded; crop_intercept = crop_intercept -50;end
crop_slope = -0.3520;

btn = uibutton(fig,"state","Text","Confirm","Position",[400,50,50,50]);
slope_slider = uislider(fig ,"Value",crop_slope,"Limits",[-1,1]);
slope_slider.Position = [100,100,250,3];
uilabel(fig,"Text","Slope","Position",[30,75,50,50])
intercept_slider = uislider(fig,"Value",crop_intercept,"Limits",[10,300]);
intercept_slider.Position = [100,50,250,3];
uilabel(fig,"Text","Intercept","Position",[30,20,50,50])

f = @(x) intercept_slider.Value-slope_slider.Value*x;
line_h = plot( ax1,ax1.YLim,f(ax1.YLim),'r','LineWidth',3);

slope_slider.ValueChangedFcn = @(sld, event) updateSlope(sld,intercept_slider,line_h);
intercept_slider.ValueChangedFcn = @(sld, event) updateIntercept(sld,slope_slider,line_h);
while ~btn.Value
    pause(.1)
end
crop_intercept = intercept_slider.Value;
crop_slope = slope_slider.Value;

%%
% fill in template scalp
template_image_filled = [];
for i = 1:im_dim(3)
    template_image_filled(:,:,i) = bwconvhull(squeeze(template_image(:,:,i)));
end
if crop
    % Define line to crop
    for i = 1:im_dim(2)
        template_image_filled(:,i,1:(im_dim(3)-ceil(-crop_slope*i+crop_intercept)))=0;
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
%% save filled template image
template_mri = mri;
template_mri.anatomy = double(template_image_filled);
template_filled_nii = [einscanpath filename(1:7) '_template_img'];
ft_write_mri([template_filled_nii,'.nii'],template_mri,'dataformat','nifti');

%% FLIRT
s = split(filename,'.');
sub = s{1};
disp("Running FSL FLIRT!");
% first command - using FLIRT to warp template scalp to einscan scalp
str1 = ['flirt -in ',template_filled_nii,'.nii -ref ',einscan_filled_nii,' -out ',template_filled_nii,'_crg.nii -omat ',template_filled_nii,'_crg.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 -interp trilinear'];
% second command - applying transformation from first command to the padded
% template mri
if padded
    str2 = ['flirt -in ',padded_mri,' -applyxfm -init ',template_filled_nii,'_crg.mat -out ',sub,'_template_file_crg.nii -paddingsize 0.0 -interp trilinear -ref ',padded_mri];
else
    str2 = ['flirt -in ',mri_file_path,' -applyxfm -init ',template_filled_nii,'_crg.mat -out ',sub '_template_file_crg.nii -paddingsize 0.0 -interp trilinear -ref ',mri_file_path];
end
stdout = system(str1); 
stdout2 = system(str2);
if stdout ~= 0 || stdout2 ~= 0
    error('Error using flirt')
end
disp("Completed FLIRT without errors")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
function updateSlope(sld,intercept_sld,line_h)
f_ = @(x) intercept_sld.Value-sld.Value*x;
line_h.YData = f_(line_h.XData);
drawnow
end
function updateIntercept(sld,slope_sld,line_h)
f_ = @(x) sld.Value-slope_sld.Value*x;
line_h.YData = f_(line_h.XData);
drawnow
end
