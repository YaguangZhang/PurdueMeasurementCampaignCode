% SIMULATESLIDINGCORRELATORCHANNELSOUNDER Simulate a sliding correlator
% channel sounder according to:
%
%   R. J. Pirkl and G. D. Durgin, "Optimal Sliding Correlator Channel
%   Sounder Design," in IEEE Transactions on Wireless Communications, vol.
%   7, no. 9, pp. 3488-3497, September 2008.
%
%   doi: 10.1109/TWC.2008.070278
%
% Yaguang Zhang, Purdue, 03/04/2019

close all; clc; clear;

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

%% Parameters

% Configure paths.
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_PROJECT_FOLDER, ...
    'ProcessingResults', '0_ChannelSounderSimulations');

% For pseudonoise (PN) signal x(t).
N = 2047;               % PN sequence length.
R_C_TX = 400e6;         % An integer chip rate at the TX side in Hz.
R_C_RX = 399.95e6;      % An integer chip rate at the RX side in Hz.
V_0 = 1;                % Height of the bipolar PN signal in volt.

% For signal simulation.
F_SIM = 20*R_C_TX;       % Simulation sample rate in Hz.

% For simulating the tapped delay line model.
TDLTotalWidthInS = 150.*10^(-9);
TDLTapWidthInS = TDLTotalWidthInS/10;
TDLDecayExp = -3;                   % Controls the extenuation rate.

% % For the millimeter wave (mm-wave) signal. 
%  F_C = 28e9;  % Carrier Frequency in Hz.

% Seed for simulation.
SIM_SEED = 999;

%% Configurations

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_PLOTS, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_PLOTS);
else
    % Clear all figures.
    previousPngFils ...
        = rdir(ABS_PATH_TO_SAVE_PLOTS, 'regexp(name, ''\d\_.+\.png'')');
    arrayfun(@(f) delete(f.name), previousPngFils);
end

% Periods.
T_C_TX = 1./R_C_TX;     % Chip period at the TX side in s.
T_C_RX = 1./R_C_RX;     % Chip period at the RX side in s.
T_SIM = 1./F_SIM;       % Simulation time step size in s.

% Slide factor.
gamma = R_C_TX/(R_C_TX-R_C_RX);

% Figure counter.
figCnt = 0;

% Some custom colors.
lightGrey = 0.9.*ones(1,3);
grey = 0.8.*ones(1,3);
darkGrey = 0.7.*ones(1,3);

% Set randam number generator state.
rng(SIM_SEED);

%% PN Sequence

% Generate the PN sequence a(i) for i = 1 to N, where a(i) is 0 or 1.
baseVal = 2;
try
    powerVal = log(N+1)/log(2);
    a = (mseq(baseVal,powerVal)+1)/2;
catch
    error(['N = ', num2str(N), ...
        ' is not valid for generating m-sequencies!']);
end

% Plot.
figCnt = figCnt+1; curFigName = 'pnSequence';
hPnSeq = figure('name', curFigName);
plotPwmPulse(a);
xlabel('Chip Index n');
ylabel('PN Sequence a(n)');

saveas(hPnSeq, fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [num2str(figCnt), '_', curFigName, '.png']));

%% PN Signal

% Convert the PN sequence a to the PN signal x. We will only observe the PN
% signal at the TX side here. Similar results hold for the PN signal at the
% RX side.
x_t_tx = @(t) V_0.*(2.* a(floor(mod(t./T_C_TX,N))+1) -1);

% One period segment of the PN signal x.
xSegTimeLengthTx = N.*T_C_TX;
xSegTimePts = 0:T_SIM:xSegTimeLengthTx;
if xSegTimePts(end) == xSegTimeLengthTx
    xSegTimePts(end) = [];
end
xSegTx = x_t_tx(xSegTimePts);

% Plot.
figCnt = figCnt+1; curFigName = 'pnSignal';
hPnSigSeg = figure('name', curFigName); hold on;
plot(xSegTimePts, xSegTx, '-', 'Color', grey);
plot(xSegTimePts, xSegTx, '.b');
xlabel('Time (s)');
ylabel('PN Signal Amplitude at TX (V)');
axis tight;
curAxis = axis;
deltaYToExtendEachSide = (curAxis(4)-curAxis(3))/10;
axis([curAxis(1:2) ...
    curAxis(3)-deltaYToExtendEachSide curAxis(4)+deltaYToExtendEachSide]);
grid minor;

saveas(hPnSigSeg, fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [num2str(figCnt), '_', curFigName, '.png']));

%% PN Signal Autocorrelation

% The autocorrelation of the PN signal can be gotten by the circular
% autocorrelation (implemented here by circular convolution) of a period
% segment of the PN signal.
numPtInXSegTx = length(xSegTx);
xTxAutoCorr = fftshift(cconv(xSegTx, xSegTx(end:-1:1), numPtInXSegTx));
halfNumPtInXSegTx = (numPtInXSegTx-1)/2;
xSegTimePtsShifted = (-halfNumPtInXSegTx:halfNumPtInXSegTx) ...
    .*T_SIM;

% Plot.
figCnt = figCnt+1; curFigName = 'pnSignalAutocorrelation';
hPnSigAutoCorr = figure('name', curFigName); hold on;
plot(xSegTimePtsShifted, xTxAutoCorr, '-', 'Color', grey);
plot(xSegTimePtsShifted, xTxAutoCorr, '.b');
xlabel('Time (s)');
ylabel('PN Signal Autocorrelation (V^2 \times s)');
axis tight;
curAxis = axis;
deltaYToExtendEachSide = (curAxis(4)-curAxis(3))/10;
axis([curAxis(1:2) ...
    curAxis(3)-deltaYToExtendEachSide curAxis(4)+deltaYToExtendEachSide]);
grid minor;

saveas(hPnSigAutoCorr, fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [num2str(figCnt), '_', curFigName, '.png']));

%% PN Signal Spectrum

% We will estimate the spectrum of the PN signal by discrete Fourier
% transfer.
xSegTxSpectrum = fftshift(fft(xSegTx));
fftLength = length(xSegTxSpectrum);
xSegTxSpectrumMag = abs(xSegTxSpectrum);
xSegTxSpectrumPha = unwrap(angle(xSegTxSpectrum));
halfFftLength = (fftLength-1)/2;
xSegTxSpectrumFPts = (-halfFftLength:halfFftLength)*(F_SIM/fftLength);

% Plot.
figCnt = figCnt+1; curFigName = 'pnSignalSpectrum';
hPnSigSpectrum = figure('name', curFigName);

subplot(2,1,1)
plot(xSegTxSpectrumFPts, xSegTxSpectrumMag)
ylabel('Magnitude (V)');
xlabel('Frequency (Hz)');
grid minor;

subplot(2,1,2)
plot(xSegTxSpectrumFPts, xSegTxSpectrumPha.*180./pi)
ylabel('Phase (Degree)');
xlabel('Frequency (Hz)');
grid minor;

saveas(hPnSigSpectrum, fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [num2str(figCnt), '_', curFigName, '.png']));

%% Sliding Correlator Spectrum

disp('Simulating the sliding correlator ...');

% For DFT:
%   Multiplication in the time domain
%       => Circular convolution in the frequency domain.
% We will start by getting the spectrum of the PN signal at the RX,
% following the procedure we used above. However, special attention needs
% to be paid to make sure the simulated PN signals at the TX and the RX
% sides have enough cycles for them to match in time while still represent
% the orignal signals correctly.

r_c_gcd = gcd(R_C_TX, R_C_RX);
numSeqTx = R_C_TX/r_c_gcd;
xSimTimeLength = N.*T_C_TX.*numSeqTx;

xSimTimePts = 0:T_SIM:xSimTimeLength;

xPnSigTx = x_t_tx(xSimTimePts);
x_t_rx = @(t) V_0.*(2.* a(floor(mod(t./T_C_RX,N))+1) -1);
xPnSigRx = x_t_rx(xSimTimePts);

% We will get the spectrum of the multiplied PN signals in two ways:
%   1. FFT the time domain signal;
%    2. Circularly convolve the spectrums.

disp('    The first method: FFT the time domain signal.')
tic;
%  % Because the computation is very time-consuming, we will cache them for
%  future usage.
% fullPathToCacheSimResults = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
%     'SimResultCache.mat');
% if exist(fullPathToCacheSimResults, 'file')
%     load(fullPathToCacheSimResults, ySpectrumFromTimeDomain);
% else
y = xPnSigTx.*xPnSigRx;
ySpectrumFromTimeDomain = fftshift(fft(y));
% save(fullPathToCacheSimResults, '-v7.3', ...
%     'y', 'ySpectrumFromTimeDomain');
% end

fftLength = length(ySpectrumFromTimeDomain);

ySpectrumFromTimeDomainMag = abs(ySpectrumFromTimeDomain);
ySpectrumFromTimeDomainPha = unwrap(angle(ySpectrumFromTimeDomain));

halfFftLength = (fftLength-1)/2;
ySpectrumFromTimeDomainFPts ...
    = (-halfFftLength:halfFftLength)*(F_SIM/fftLength);

% We will generate two figures for the spectrum, one overview for all
% frequency range and one zoomed in for details.
decNumForOverview = 100;
zoomedInFreqRange = [-5, 5].*(R_C_TX-R_C_RX);

boolsPtsToShowZoomed ...
    = (ySpectrumFromTimeDomainFPts>=zoomedInFreqRange(1)) ...
    & (ySpectrumFromTimeDomainFPts<=zoomedInFreqRange(2));

% Plot.
figCnt = figCnt+1;
curFigName = 'multipliedPnSignalSpectrumViaDftTimeDomainSigOverview';
hMultipliedPnSigSpectrumFromTimeOverview = figure('name', curFigName);

subplot(2,1,1)
plot(ySpectrumFromTimeDomainFPts(1:decNumForOverview:end), ...
    ySpectrumFromTimeDomainMag(1:decNumForOverview:end));
ylabel('Magnitude (V \times s)');
xlabel('Frequency (Hz)');
grid minor;

subplot(2,1,2)
plot(ySpectrumFromTimeDomainFPts(1:decNumForOverview:end), ...
    ySpectrumFromTimeDomainPha(1:decNumForOverview:end).*180./pi);
ylabel('Phase (Degree)');
xlabel('Frequency (Hz)');
grid minor;

saveas(hMultipliedPnSigSpectrumFromTimeOverview, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [num2str(figCnt), '_', curFigName, '.png']));

% Plot.
figCnt = figCnt+1;
curFigName = 'multipliedPnSignalSpectrumViaDftTimeDomainSigZoomed';
hMultipliedPnSigSpectrumFromTimeZoomed = figure('name', curFigName);

subplot(2,1,1)
plot(ySpectrumFromTimeDomainFPts(boolsPtsToShowZoomed), ...
    ySpectrumFromTimeDomainMag(boolsPtsToShowZoomed));
ylabel('Magnitude (V \times s)');
xlabel('Frequency (Hz)');
grid minor;

subplot(2,1,2)
plot(ySpectrumFromTimeDomainFPts(boolsPtsToShowZoomed), ...
    ySpectrumFromTimeDomainPha(boolsPtsToShowZoomed).*180./pi);
ylabel('Phase (Degree)');
xlabel('Frequency (Hz)');
grid minor;

saveas(hMultipliedPnSigSpectrumFromTimeZoomed, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [num2str(figCnt), '_', curFigName, '.png']));

toc;

% Clear huge variables to free memory.
close all;
clearvars ySpectrumFromTimeDomainMag ySpectrumFromTimeDomainPha;

disp(' ');
disp('    The second method: Circularly convolve the spectrums.');
tic;

xPnSigTxSpectrum = fft(xPnSigTx);
xPnSigRxSpectrum = fft(xPnSigRx);

ySpecturmFromFreqDomain ...
    = cconv(xPnSigTxSpectrum, xPnSigRxSpectrum, fftLength);
ySpecturmFromFreqDomain = fftshift(ySpecturmFromFreqDomain);

assert(fftLength == length(ySpecturmFromFreqDomain), ...
    'Spectrums computed using different methods should have the same length!');

% Clear huge variables to free memory.
clearvars xPnSigTxSpectrum xPnSigRxSpectrum;
ySpectrumFromFreqDomainMag = abs(ySpecturmFromFreqDomain);
ySpectrumFromFreqDomainPha = unwrap(angle(ySpecturmFromFreqDomain));

% Plot.
figCnt = figCnt+1;
curFigName = 'multipliedPnSignalSpectrumViaCirConvSpecsOverview';
hMultipliedPnSigSpectrumFromFreqOverview = figure('name', curFigName);

subplot(2,1,1)
plot(ySpectrumFromTimeDomainFPts(1:decNumForOverview:end), ...
    ySpectrumFromFreqDomainMag(1:decNumForOverview:end));
ylabel('Magnitude (V \times s)');
xlabel('Frequency (Hz)');
grid minor;

subplot(2,1,2)
plot(ySpectrumFromTimeDomainFPts(1:decNumForOverview:end), ...
    ySpectrumFromFreqDomainPha(1:decNumForOverview:end).*180./pi);
ylabel('Phase (Degree)');
xlabel('Frequency (Hz)');
grid minor;

saveas(hMultipliedPnSigSpectrumFromFreqOverview, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [num2str(figCnt), '_', curFigName, '.png']));

% Plot.
figCnt = figCnt+1;
curFigName = 'multipliedPnSignalSpectrumViaCirConvSpecsZoomed';
hMultipliedPnSigSpectrumFromFreqZoomed = figure('name', curFigName);

subplot(2,1,1)
plot(ySpectrumFromTimeDomainFPts(boolsPtsToShowZoomed), ...
    ySpectrumFromFreqDomainMag(boolsPtsToShowZoomed));
ylabel('Magnitude (V \times s)');
xlabel('Frequency (Hz)');
grid minor;

subplot(2,1,2)
plot(ySpectrumFromTimeDomainFPts(boolsPtsToShowZoomed), ...
    ySpectrumFromFreqDomainPha(boolsPtsToShowZoomed).*180./pi);
ylabel('Phase (Degree)');
xlabel('Frequency (Hz)');
grid minor;

saveas(hMultipliedPnSigSpectrumFromFreqZoomed, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [num2str(figCnt), '_', curFigName, '.png']));

% Clear unnecessary variables to free up memory.
clearvars ySpectrumFromFreqDomainMag ySpectrumFromFreqDomainPha;

toc;

disp('Done!');

%% Low-Pass Filter (LPF)

disp('Simulating filtering the multiplication result ...');

tic;

% We will reuse the result before to save memory.
%    y = xPnSigTx.*xPnSigRx; ySpectrumFromTimeDomain = fft(y);
ySpectrumFromTimeDomain = ifftshift(ySpectrumFromTimeDomain);

% Use the optimal LPF.
LPFFreqRange = [-1 1].*(R_C_TX-R_C_RX);

% Note that the FFT results are stored for frequency range [0, F_SIM) Hz.
fftLength = length(ySpectrumFromTimeDomain);
halfFftLength = (fftLength-1)/2;
ySpectrumFromTimeDomainFPts ...
    = (-halfFftLength:halfFftLength)*(F_SIM/fftLength);
boolsLPEliminated = ifftshift( ...
    (ySpectrumFromTimeDomainFPts<LPFFreqRange(1)) ...
    | (ySpectrumFromTimeDomainFPts>LPFFreqRange(2)));
ySpectrumFromTimeDomainLowPassed = ySpectrumFromTimeDomain;
ySpectrumFromTimeDomainLowPassed(boolsLPEliminated) = 0;
ySpectrumFromTimeDomainHighPassed = ySpectrumFromTimeDomain;
ySpectrumFromTimeDomainHighPassed(~boolsLPEliminated) = 0;

% For displaying the spectrum of the passed signal.
ySpectrumFromTimeDomainLowPassedShiftedToShow ...
    = fftshift(ySpectrumFromTimeDomainLowPassed);
ySpectrumFromTimeDomainLowPassedShiftedToShow ...
    = ySpectrumFromTimeDomainLowPassedShiftedToShow(boolsPtsToShowZoomed);
ySpectrumFromTimeDomainLowPassedMagToShow ...
    = abs(ySpectrumFromTimeDomainLowPassedShiftedToShow);
ySpectrumFromTimeDomainLowPassedPhaToShow ...
    = unwrap(angle(ySpectrumFromTimeDomainLowPassedShiftedToShow));

yLowPassed = ifft(ySpectrumFromTimeDomainLowPassed);
yHighPassed = y-yLowPassed;

% We will shift and down sample the result in the plot.
decNumForOverview = 100;
numSamps = length(y);
halfNumSamps = (numSamps-1)/2;

yToShow = fftshift(y);
yLowPassedToShow = fftshift(yLowPassed);
yHighPassedToShow = fftshift(yHighPassed);
timePts = (-halfNumSamps:halfNumSamps)'./F_SIM;

yToShow = yToShow(1:decNumForOverview:end);
yLowPassedToShow = yLowPassedToShow(1:decNumForOverview:end);
yHighPassedToShow = yHighPassedToShow(1:decNumForOverview:end);
timePtsToShow = timePts(1:decNumForOverview:end);

% We will show three periods for better results.
numPeriodsToShow = 3;
yToShow = repmat(yToShow, numPeriodsToShow, 1);
yLowPassedToShow = repmat(yLowPassedToShow, numPeriodsToShow, 1);
yHighPassedToShow = repmat(yHighPassedToShow, numPeriodsToShow, 1);
timeShift = numSamps./F_SIM;

timeRangeToShow = [timePtsToShow(1)-0.6*timeShift, ...
    timePtsToShow(end)+0.6*timeShift];
timePtsToShow = [timePtsToShow-timeShift; ...
    timePtsToShow; ...
    timePtsToShow+timeShift];

% Plot.
figCnt = figCnt+1;
curFigName = 'multipliedPnSignalLowPassedVsDejected';
hMultipliedPnSignalLowPassed = figure('name', curFigName);
hold on;

hHP = plot(timePtsToShow, yHighPassedToShow, '-', 'Color', darkGrey);
hLP = plot(timePtsToShow, yLowPassedToShow, '-.b');
ylabel('Magnitude (V^2)');
xlabel('Time (s)');
legend([hLP, hHP], 'Low Passed', 'Rejected');
axis([timeRangeToShow, -2 2]);
grid minor;

saveas(hMultipliedPnSignalLowPassed, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [num2str(figCnt), '_', curFigName, '.png']));

% Plot.
figCnt = figCnt+1;
curFigName = 'multipliedPnSignalLowPassedVsDejectedZoomedIn';
hMultipliedPnSignalLowPassedZoomed = figure('name', curFigName);
hold on;

hHP = plot(timePtsToShow, yHighPassedToShow, '-', 'Color', darkGrey);
hLP = plot(timePtsToShow, yLowPassedToShow, '-.b');
ylabel('Magnitude (V^2)');
xlabel('Time (s)');
legend([hLP, hHP], 'Low Passed', 'Rejected');
axis([timeRangeToShow.*0.005, -2 2]);
grid minor;

saveas(hMultipliedPnSignalLowPassedZoomed, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [num2str(figCnt), '_', curFigName, '.png']));

% Plot.
figCnt = figCnt+1;
curFigName = 'multipliedPnSignalLowPassedSpectrum';
hMultipliedPnSignalLowPassedSpectrum = figure('name', curFigName);

freqRangeToShow = LPFFreqRange.*1.5;

subplot(2,1,1)
plot(ySpectrumFromTimeDomainFPts(boolsPtsToShowZoomed), ...
    ySpectrumFromTimeDomainLowPassedMagToShow);
ylabel('Magnitude (V \times s)');
xlabel('Frequency (Hz)');
xlim(freqRangeToShow);
grid minor;

subplot(2,1,2)
plot(ySpectrumFromTimeDomainFPts(boolsPtsToShowZoomed), ...
    ySpectrumFromTimeDomainLowPassedPhaToShow.*180./pi);
ylabel('Phase (Degree)');
xlabel('Frequency (Hz)');
xlim(freqRangeToShow);
grid minor;

saveas(hMultipliedPnSignalLowPassedSpectrum, ...
    fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    [num2str(figCnt), '_', curFigName, '.png']));

% Clear unnecessary variables to free up memory.
clearvars y yHighPassed yHighPassedToShow yLowPassed yLowPassedToShow ...
    ySpectrumFromTimeDomainLowPassed ySpectrumFromTimeDomainHighPassed  ...
    ySpectrumFromTimeDomainLowPassedMag ...
    ySpectrumFromTimeDomainLowPassedPha ...
    ySpectrumFromTimeDomainLowPassedShifted;

toc;

disp('Done!');

%% Transmission through the Tapped Delay Line Model
% We will simulate multiple implementations of the tapped delay line model.

disp('Simulating the sliding correlator with the tapped delay line model ...');

numTDLSims = 5;
maxNumTaps = 5;
for dixTDLSim = 1:numTDLSims
    disp(['    Simulation ', ...
        num2str(dixTDLSim), '/', num2str(numTDLSims), '...']);
    tic;
    
    % Simulate the channel.
    curNumTaps = 1+floor(rand(1).*maxNumTaps);
    [curTDLImpulseResp, curTDLImpulseRespTimePts] ...
        = genTDLImpulseResponse(curNumTaps, F_SIM, ...
        TDLTotalWidthInS, TDLTapWidthInS, TDLDecayExp);
    
    % We will FFT the time domain signal because it is faster.
    yTDL = cconv(xPnSigTx, curTDLImpulseResp, length(xPnSigTx)).*xPnSigRx;
    assert(length(yTDL) == numSamps, 'Variable numSamps is outdated!');
    
    % LPF the RX signal and covert it back to time domain.
    yTDLSpectrum = fft(yTDL);
    yTDLSpectrum(boolsLPEliminated) = 0;
    yTDL = ifft(yTDLSpectrum);
    
    % Down sample signals for plotting.
    yTDLToShow = fftshift(yTDL(1:decNumForOverview:end));
    yTDLTimePtsToShow = timePts(1:decNumForOverview:end)/gamma;
    
    % Align the RX signal with the channel impulse response according to
    % their maximum value location.
    [yTDLMax, yTDLArgMaxIdx] = max(yTDLToShow);
    [channelMax, channelArgMaxIdx] = max(curTDLImpulseResp);
    
    peakTime = curTDLImpulseRespTimePts(channelArgMaxIdx);
    yTDLTimePtsToShow = yTDLTimePtsToShow ...
        + peakTime - yTDLTimePtsToShow(yTDLArgMaxIdx);
    
    % Plot.
    figCnt = figCnt+1;
    curFigName = 'dilutedOutputVsChannel';
    hDilutedOutputVsChannel = figure('name', curFigName);
    hold on;
    
    colorChannel = 'k';
    colorDilutatedOutput = 'b';
    
    hChannel = plot( ...
        curTDLImpulseRespTimePts, curTDLImpulseResp/channelMax, ...
        '--', 'Color', colorChannel);
    axC = gca; % current axes
    axC.XColor = colorChannel;
    axC.YColor = colorChannel;

    hDilutedOutput = plot(yTDLTimePtsToShow, yTDLToShow/yTDLMax, ...
        '-', 'Color', colorDilutatedOutput);
    
    ylabel('Normalized Magnitude');
    xlabel('Time (s)');
    legend([hDilutedOutput, hChannel], 'Time-Diluted Output', 'Channel');
    grid minor;
    axis tight;
    xlim([curTDLImpulseRespTimePts(1), curTDLImpulseRespTimePts(end)]);
    
    saveas(hDilutedOutputVsChannel, ...
        fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        [num2str(figCnt), '_', curFigName, '.png']));
    
    toc;
end

disp('Done!');
% EOF