%% Housekeeping
clear all
close all
clc

addpath('C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\Braille_kids\')
addpath('C:\Users\ppanr\Documents\GitHub\Braille_kids')
project_dir = 'R:\DRS-KidsOPM\Temp_warp_paper\individual\';
%project_dir =  '/rdrives/ppzlr/DRS-KidsOPM/Paediatric_OPM_Notts/'; %linux
% %% Set subject numbers to run
%good_subs = [1:13 15 16 18 21 24:26];
%
good_subs = 18;
%% Tstat from 1 mm sensorimotor mask
for sub_i = good_subs
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
    ses_i = 1;
    ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0'
    tic
    index_pinky_Tstat_func(sub,ses,project_dir)
    toc
end

% TFS from tstat peak                        
for sub_i = good_subs
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
    ses_i = 1
    ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0'
    index_pinky_peak_TFS_func(sub,ses,project_dir)
end

% Tstat from 4 mm whole brain
for sub_i = good_subs
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
    ses_i = 1
    ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0'
    index_pinky_4mm_Tstat_func(sub,ses,project_dir)
end

% Connectivity
for sub_i = good_subs
    sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
    ses_i = 1
    ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0' 
    kids_opm_test_func(sub,ses,project_dir)
end

% Lead fields
% for sub_i = good_subs
%     sub = sprintf('1%2d',sub_i);sub(sub == ' ') = '0'
%     ses_i = 1
%     ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0' 
%     get_aal_lf(sub,ses,project_dir)
% end



