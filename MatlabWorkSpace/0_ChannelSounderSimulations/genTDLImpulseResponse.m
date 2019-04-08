function [ TDLImpulseResp, TDLImpulseRespTimePts] ...
    = genTDLImpulseResponse(numTaps, sampleRate, ...
    TDLTotalWidthInS, TDLTapWidthInS, TDLDecayExp)
%GENTDLIMPULSERESPONSE Generate an impulse response for an implementation
%of the Tapped Delay Line (TDL) model.
%
% Inputs:
%   - numTaps
%     The number of taps.
%   - sampleRate
%     The sample rate for the impulse response in Hz.
%   - TDLTotalWidthInS, TDLTapWidthInS
%     The total width of the output impulse response in time.
%   - TDLTapWidthInS
%     The width of each tap in time.
%   - TDLDecayExp
%     We will use a exponential decay to control the taps' amplitudes.
%
% Outputs:
%   - TDLImpulseResp
%     The samples (organized as a column vector) for the simulated impulse
%     response of the Tapped Delay Line (TDL) model.
%   - TDLImpulseRespTimePts
%     The sampling time for elements in TDLImpulseResp.
%
% Yaguang Zhang, Purdue, 03/04/2019

assert(numTaps>=1, 'There should be at least one tap!');

% Initialize the result.
numSamps = ceil(sampleRate.*TDLTotalWidthInS);
% We will first organize the samples for each tap individually as rows of a
% matrix.
TDLImpulseRespMat = zeros(numTaps, numSamps);

% We will control the tap's amplitude according to its location.
getTapAmplitude = @(x) exp(TDLDecayExp.*x);
% The first tap always starts at the very beginning.
tapStartLocNormed = [0; rand(numTaps-1, 1)];
tapAmpsNormed = getTapAmplitude(tapStartLocNormed);

% Generate tap samples.
numSampsPerTap = ceil(sampleRate.*TDLTapWidthInS);
numSampsPerTapHalf = floor(numSampsPerTap/2);

tapStartSampleIndices ...
    = floor(tapStartLocNormed'.*(numSamps-numSampsPerTap))+1;

% Some randomness is introduced to perturb the tap amplitudes.
randomRatioForAmplitude = 0.1;
for idxTap = 1:numTaps
    curTapAmp = tapAmpsNormed(idxTap) ...
        .*(1-randomRatioForAmplitude+rand(1).*randomRatioForAmplitude.*2);
    curTapSamps = linspace(0, curTapAmp, numSampsPerTapHalf);
    TDLImpulseRespMat(idxTap, ...
        tapStartSampleIndices(idxTap)...
        :(tapStartSampleIndices(idxTap)+2.*length(curTapSamps)-2)) ...
        = [curTapSamps curTapSamps((end-1):-1:1)];
end

% For overlapping taps, we will use the largest value.
TDLImpulseResp = max(TDLImpulseRespMat, [], 1);

% For the corresponding sampling time.
TDLImpulseRespTimePts = linspace(0, TDLTotalWidthInS, numSamps);

end
% EOF