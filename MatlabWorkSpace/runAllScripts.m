% RUNALLSCRIPTS Post process the measurement dataset.
%
% This is a holder script listing all the steps available to properly
% process the measurement dataset. Please comment & uncomment commands as
% it is needed, depending on which results require updates.
%
% Yaguang Zhang, Purdue, 03/14/2019

clear; clc; close all;

% Add libs to current path and set ABS_PATH_TO_PROJECT_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(pwd)));
setPath;

%% 0_ChannelSounderSimulations
addpath(fullfile(pwd, '0_ChannelSounderSimulations'));
SimulateSlidingCorrelatorChannelSounder;

% EOF