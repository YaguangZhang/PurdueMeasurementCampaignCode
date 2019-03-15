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

close all; clc;

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
F_SIM = 10e9;           % Simulation sample rate in Hz.
T_SIM_MAX = 0.1;        % Maximum time to simulate.

% For the millimeter wave (mm-wave) signal.
F_C = 28e9;             % Carrier Frequency in Hz.

%% Configurations

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_PLOTS, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_PLOTS);
end

% Periods.
T_C_TX = 1./R_C_TX;     % Chip period at the TX side in s.
T_C_RX = 1./R_C_RX;     % Chip period at the RX side in s.
T_SIM = 1./F_SIM;       % Simulation time step size in s.

% Figure counter.
figCnt = 0;

% Some custom colors.
lightGrey = 0.9.*ones(1,3);
grey = 0.8.*ones(1,3);
darkGrey = 0.7.*ones(1,3);

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
ylabel('PN Signal Autocorrelation (V^2)');
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

% For DFT:
%   Multiplication in the time domain
%       => Circular convolution in the frequency domain.
% We will start by getting the spectrum of the PN signal at the RX,
% following the procedure we used above. However, special attention needs
% to be paid to make sure the simulated PN signals at the TX and the RX
% sides have enough cycles for them to match in time while still represent
% the orignal signals correctly.

[r_c_gcd, r_c_tx_min, r_c_rx_min] = gcd(R_C_TX, R_C_RX);

x_t_rx = @(t) V_0.*(2.* a(floor(mod(t./T_C_RX,N))+1) -1);
xSegRx = x_t_rx(xSegTimePts);
xSegRxSpectrum = fftshift(fft(xSegRx));
assert(fftLength == length(xSegRxSpectrum), ...
    'Unable to deal with spectrums of different sizes!');
xSegRxSpectrumMag = abs(xSegRxSpectrum);
xSegRxSpectrumPha = unwrap(angle(xSegRxSpectrum));
xSegRxSpectrumFPts = (-halfFftLength:halfFftLength)*(F_SIM/fftLength);

% EOF