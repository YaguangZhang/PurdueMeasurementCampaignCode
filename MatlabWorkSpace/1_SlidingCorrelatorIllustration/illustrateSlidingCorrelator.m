% ILLUSTRATESLIDINGCORRELATOR This script will generate vidoe clips to
% illustrate how the sliding correlator works.
%
% Yaguang Zhang, Purdue, 07/17/2021

close all; clc; clearvars -except FLAG_GEN_ANIME N;

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

%% Parameters

% Configure paths.
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_PROJECT_FOLDER, ...
    'ProcessingResults', '1_SlidingCorrelatorIllustration');

if ~exist('FLAG_GEN_ANIME', 'var')
    FLAG_GEN_ANIME = true;
end

% For pseudonoise (PN) signal x(t).
if ~exist('N', 'var')
    N = 7;              % PN sequence length (2^k-1).
end
R_C_TX = N;             % An integer chip rate at the TX side.
R_C_RX = N-1;           % An integer chip rate at the RX side
V_0 = 1;                % Height of the bipolar PN signal in volt.

% Seed for simulation.
SIM_SEED = 999;

%% Configurations

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_PLOTS, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_PLOTS);
end

% Some custom colors.
lightGrey = 0.9.*ones(1,3);
darkGrey = 0.7.*ones(1,3);

startLineC = 'r';
startLinePattern = '--';

normLineW = 1;
thickLineW = 2;
subPlotFontSize = 10;

% Set randam number generator state.
rng(SIM_SEED);

%% Simulation

% Generate the PN sequence a(i) for i = 1 to N, where a(i) is 0 or 1.
baseVal = 2;
try
    powerVal = log(N+1)/log(2);
    a = mseq(baseVal,powerVal);
catch
    error(['N = ', num2str(N), ...
        ' is not valid for generating m-sequencies!']);
end

% Construct the TX and RX signals for one cycle.
totalTimeSteps = R_C_TX*R_C_RX;
[sigTx, sigRx] = deal(nan(totalTimeSteps, 1));
for idxSample = 1:N
    sigTx((idxSample-1).*R_C_RX + (1:R_C_RX)) = a(idxSample);
    sigRx((idxSample-1).*R_C_TX + (1:R_C_TX)) = a(idxSample);
end

% Repeat these signals until they are aligned again.
sigTxCycleLen = length(sigTx);
sigRxCycleLen = length(sigRx);
sigTx = repmat(sigTx, [R_C_TX, 1]);
sigRx = repmat(sigRx, [R_C_RX, 1]);
assert(length(sigTx)==length(sigRx), ...
    'The extended TX and RX signals are not of the same length!')

sigTxDoubled = [sigTx; sigTx];
sigRxDoubled = [sigRx; sigRx];

% Doing the correlation with the TX signal as the reference.
simTimeLength = length(sigTx);
sliCor = sigTx.*nan;
correlationWindowSize = sigTxCycleLen;
for idxT = correlationWindowSize:(simTimeLength+correlationWindowSize-1)
    curIndices = (idxT-correlationWindowSize+1):idxT;
    sliCor(idxT) = sum( ...
        sigTxDoubled(curIndices).*sigRxDoubled(curIndices) ...
        )/correlationWindowSize;
end
sliCor(1:(correlationWindowSize-1)) ...
    = sliCor((end-correlationWindowSize+2):end);
sliCor((length(sigTx)+1):end) = [];

assert(length(sigTx)==length(sliCor), ...
    ['The TX signal and correlation results ', ...
    'are not of the same length!'])

%% Visualization

NUM_OF_TX_SIG_CYCLE_TO_SHOW = 2;

snipetToPlotOverview;
hOverviewForOneAlignment = curFig;
xlim([1, simTimeLength]);

saveas(hOverviewForOneAlignment, fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    ['AnimationForOneAlignment_N_', num2str(N), '.jpg']));

%% Animation

if FLAG_GEN_ANIME
    QUALITY = 95;
    FRAME_RATE = 60;
    VIDEO_LEN_IN_S = 10;
    
    snipetToPlotOverview;
    hAnimeForOneAlignment = curFig;
    
    totalNumOfFrames = FRAME_RATE*VIDEO_LEN_IN_S;
    simTsToGenNewFrame = floor(linspace(sigTxCycleLen, ...
        simTimeLength+sigTxCycleLen-1, ...
        totalNumOfFrames));
    
    v = VideoWriter(fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        ['AnimationForOneAlignment_N_', num2str(N), '.avi']), ...
        'Motion JPEG AVI');
    v.Quality = QUALITY;
    v.FrameRate = FRAME_RATE;
    open(v);
    
    for curT = simTsToGenNewFrame
        xlim([curT-sigTxCycleLen+1, curT]);
        frame = getframe(hAnimeForOneAlignment);
        writeVideo(v,frame);
    end
    close(v);
end

% EOF