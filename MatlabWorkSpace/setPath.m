% SETPATH Add lib folders into Matlab path.
%
% Yaguang Zhang, Purdue, 03/14/2019

cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
addpath(genpath(fullfile(pwd, 'lib')));

% The absolute path to the shared folder holding the data and code for the
% Purdue mm-wave data. Please make sure it is correct for the machine which
% will run this script.
%  - On (quite powerful) Windows Artsy:
absPathWinArtsy = 'D:\One Drive - Purdue\OneDrive - purdue.edu\EARS\2019_Purdue Measurement Campaign';
%  - Local copy on Windows Dell:
absPathWinDell = 'C:\Users\Zyglabs\OneDrive - purdue.edu\EARS\2019_Purdue Measurement Campaign';
unknownComputerErrorMsg = ...
    ['Compute not recognized... \n', ...
    '    Please update setPath.m for your machine. '];
unknownComputerErrorId = 'setPath:computerNotKnown';
switch getenv('computername')
    case 'ZYGLABS-DELL'
        % ZYG's Dell laptop.
        ABS_PATH_TO_PROJECT_FOLDER = absPathWinDell;
    case 'ARTSY'
        % ZYG's lab desktop.
        ABS_PATH_TO_PROJECT_FOLDER = absPathWinArtsy;
    otherwise
        error(unknownComputerErrorId, unknownComputerErrorMsg);
end

%EOF