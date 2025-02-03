function get_aal_lf(sub,ses,project_dir)
restoredefaultpath
cleaning_only = 0;
close all
clc
% project_dir = 'P:\Paediatric_OPM_Notts\';
% sub = '006';
% ses = '001';
run = 'run-001'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath([project_dir,filesep,'fieldtrip-20220906'])
addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
addpath('C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\Braille_kids\Beamformer')
addpath('D:\GitHub\BrailleKids\Beamformer\')
script_dir = mfilename('fullpath');fname = mfilename;script_dir = script_dir(1:end-length(fname));
addpath(script_dir)
addpath([script_dir,'Beamformer',filesep,''])
ft_defaults;
datadir = [project_dir,'Data',filesep,'BIDS',filesep];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_type = '_task-braille';
filename = ['sub-',sub,'_ses-',ses,exp_type,'_',run];
path_main = [datadir,'sub-',sub,filesep,'ses-',ses,filesep];

path_ICA = [datadir,'derivatives',filesep,'ICA',filesep,'sub-',sub,filesep];
if ~exist(path_ICA,'dir'); mkdir(path_ICA);end
files_ICA = [path_ICA,filename];

path_cleaning = [datadir,'derivatives',filesep,'cleaning',filesep,'sub-',sub,filesep];
if ~exist(path_cleaning,'dir'); mkdir(path_cleaning);end
files_cleaning = [path_cleaning,filename];

path_VEs = [datadir,'derivatives',filesep,'VEs',filesep,'sub-',sub,filesep];
if ~exist(path_VEs,'dir'); mkdir(path_VEs);end
files_VEs = [filename,'_VE'];

path_AEC = [datadir,'derivatives',filesep,'AEC',filesep,'sub-',sub,filesep];
if ~exist(path_AEC,'dir'); mkdir(path_AEC);end
files_AEC = [filename,'_AEC'];

path_meg_data = [path_main,'meg',filesep];

path_meshes = [datadir,'derivatives',filesep,'sourcespace',filesep,'sub-',sub,filesep];
files_meshes = ['sub-',sub,'_meshes.mat'];
files_AAL_centroids = ['sub-',sub,'_AAL_centroids.nii.gz'];
files_AAL_regions = ['sub-',sub,'_AAL_regions.nii.gz'];

path_mri = [path_main,'anat',filesep];
files_mri = ['sub-',sub,'_anat.nii'];
S.mri_file = [path_mri,files_mri]; 

path_helmet = [datadir,'derivatives',filesep,'helmet',filesep,'sub-',sub,filesep];
files_helmet_info = dir([path_helmet,'*.mat']);files_helmet_info=files_helmet_info.name;

files_voxlox = ['sub-',sub,'_AAL_locs.mat'];
files_channels = [filename,'_channels.tsv'];

files_events = [filename,'_events.tsv'];


%read meg data
cd(path_meg_data)
read_info = readlines([filename,'_meg_read_info.txt']);
Size = strsplit(read_info(1));Size = [str2num(Size(2)),str2num(Size(4))];

Precision = strsplit(read_info(2));Precision = Precision(2);

Ordering = strsplit(read_info(3));Ordering = Ordering(2);

FileID = fopen([path_meg_data,filename,'_meg.dat'],'r');

data=fread(FileID,Size,lower(Precision),Ordering)';
fclose(FileID);

% fs and other info
fID = fopen([path_meg_data,filename,'_meg.json']);
raw = fread(fID,inf);
json_info = jsondecode(char(raw'));
fs = json_info.SamplingFrequency;

% trigger info
event_table = readtable([path_meg_data,filename,'_events.tsv'],'FileType','text','delimiter','\t');
all_start_samps = [round(event_table(startsWith(event_table.type(:),'Start_index'),:).sample);...
    round(event_table(startsWith(event_table.type(:),'Start_pinky'),:).sample)];

%helmet info
load([path_helmet,files_helmet_info])

%% Preproc
% Mean correct
data = data - mean(data,1);

% Notch filter
for harms = [50,100,150]
    Wo = harms/(fs/2);  BW = Wo/35;
    [b,a] = iirnotch(Wo,BW);
    disp('Applying Notch filter')
    data = filter(b,a,data,[],1);
end
Nchans = size(data,2);

% bandpass filter for viewing
disp('Applying 1-150 Hz bandpass filter')

hp = 1;
lp = 150;
[b,a] = butter(4,2*[hp lp]/fs);
data_f = [filtfilt(b,a,data)]';

% Get rid of bad channels
%% Get rid of bad channels

if ~exist([path_meg_data,files_channels(1:end-4),'_proc.tsv'],'file')
%     error("Use 'get_good_channels.m for this dataset first")
    ch_table = readtable([path_meg_data,files_channels],'FileType','text','Delimiter','tab');
    ch_table.isx = endsWith(ch_table.name,'X');
    ch_table.isy = endsWith(ch_table.name,'Y');
    ch_table.isz = endsWith(ch_table.name,'Z');
    ch_table.slot_no = zeros(height(ch_table),1);
    % sanity check
    if sum(sum([ch_table.isx,ch_table.isy,ch_table.isz],2)) ~= height(ch_table)
        error('Channel orientation [x,y,z] labels might be wrong!')
    end
    [ch_table] = Bad_Channels(data_f',ch_table,fs);
    writetable(ch_table,[path_meg_data,files_channels(1:end-4),'_proc.tsv'],...
        'WriteRowNames',true,'Delimiter','tab','FileType','text')
else
    ch_table = readtable([path_meg_data,files_channels(1:end-4),'_proc.tsv'],...
        'Delimiter','tab','FileType','text');
end

%% getting slot numbers for later. Change ICA comp visualiser?
precision_limit = 1e-6;
pos = [ch_table.Px,ch_table.Py,ch_table.Pz];
for sl_i = 1:size(Helmet_info.lay.pos,1)
    detected_ind = find(sqrt(sum((repmat(Helmet_info.sens_pos(sl_i,:),height(ch_table),1) - pos/100).^2,2)) < precision_limit)
    if ~isempty(detected_ind)
        ch_table.slot_no(detected_ind) = sl_i;
    end
end
%%
% remove from data matrix
disp("Removing bad channels")
bad_chans_data = find(startsWith(ch_table.status,'bad'));

ch_table(bad_chans_data,:) = [];
data_f(bad_chans_data,:) = [];

%% sensor info
S.sensor_info.pos = [ch_table.Px,ch_table.Py,ch_table.Pz]./100;
S.sensor_info.ors = [ch_table.Ox,ch_table.Oy,ch_table.Oz];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

xfm = load(sprintf('sub-%s_ses-%s%s_%s_sens2template_transform.txt',sub,ses,exp_type,run));
xfm = reshape(xfm,4,4);
xfm_rot = xfm(1:3,1:3);
xfm_translation = xfm(1:3,4)';
S.sensor_info.pos = (xfm_rot*S.sensor_info.pos' + xfm_translation')';
S.sensor_info.ors = (xfm_rot*S.sensor_info.ors')';

%% Source positions

% load AAL locations 
AAL_locs_mat_file = [path_meshes,files_AAL_centroids(1:end-7),'.mat'];
if ~exist(AAL_locs_mat_file,'file')
    AAL_regions = ft_read_mri([path_meshes,files_AAL_regions]);
    AAL_regions = ft_convert_units(AAL_regions,'m');
    [sourcepos_vox] = get_AAL_coords(AAL_regions,S)
    sourcepos = ft_warp_apply(AAL_regions.transform,sourcepos_vox);
    save(AAL_locs_mat_file,'sourcepos')
else
    load(AAL_locs_mat_file,'sourcepos')
end
[bf_outs_shell] = run_beamformer('shell',sourcepos,S,0,[],0);
lead_fields_shell_xyz = bf_outs_shell.LF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% convert orientation of sources to polar
% Load meshes

load([path_meshes,files_meshes],'meshes');

meshes = ft_convert_units(meshes,'m');

X = meshes.pnt; Origin = mean(X,1);
Ndips = length(sourcepos);
for n = 1:Ndips
    thispos = sourcepos(n,:);
    [phi,theta1,~] = cart2sph(thispos(1) - Origin(1),thispos(2) - Origin(2) ,thispos(3) - Origin(3));
    theta = pi/2 - theta1;
    Src_Or_theta(n,:) = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)];
    Src_Or_phi(n,:) = [-sin(phi) cos(phi) 0];
    Lead_fields(:,1,n) = Src_Or_theta(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_theta(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_theta(n,3)*lead_fields_shell_xyz(:,3,n);
    Lead_fields(:,2,n) = Src_Or_phi(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_phi(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_phi(n,3)*lead_fields_shell_xyz(:,3,n);
end
save([path_AEC filename '_lead_fields.mat'],'Lead_fields','sourcepos','S');
figure;
plot3(squeeze(sourcepos(:,1)),squeeze(sourcepos(:,2)),squeeze(sourcepos(:,3)),'o')
hold on
plot3(S.sensor_info.pos(:,1),S.sensor_info.pos(:,2),S.sensor_info.pos(:,3),'kx','MarkerSize',10)
end
