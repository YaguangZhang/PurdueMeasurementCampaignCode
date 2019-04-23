function [ TDLImpulseResp, TDLImpulseRespTimePts] ...
    = genTruncatedTrainImpulseResponse(numTaps, sampleRate, ...
    TDLTapWidthInS, TDLTapDelayInS)
%GENTRUNCATEDTRAINIMPULSERESPONSE Generate a truncated train (the Tapped
%Delay Line model with the exact same shape for each tap and the same delay
%between taps) impulse response.
%
% Inputs:
%   - numTaps
%     The number of taps.
%   - sampleRate
%     The sample rate for the impulse response in Hz.
%   - TDLTapWidthInS
%     The width of each tap in time.
%   - TDLTapDelayInS
%     The delay between adjacent taps.
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
assert(TDLTapWidthInS>0, 'Tap width should be positive!');

% The total width of the output impulse response in time.
TDLTotalWidthInS = TDLTapWidthInS + TDLTapDelayInS*(numTaps-1);

% Initialize the result.
numSamps = ceil(sampleRate.*TDLTotalWidthInS);
% We will first organize the samples for each tap individually as rows of a
% matrix.
TDLImpulseRespMat = zeros(numTaps, numSamps);

% The first tap always starts at the very beginning.
tapStartLocNormed = linspace(0, 1, numTaps);
tapAmpsNormed = ones(numTaps, 1);

% Generate tap samples.
numSampsPerTap = ceil(sampleRate.*TDLTapWidthInS);
numSampsPerTapHalf = floor(numSampsPerTap/2);

tapStartSampleIndices ...
    = floor(tapStartLocNormed'.*(numSamps-numSampsPerTap))+1;

% Unlike the original TDL model, no randomness is introduced to perturb the
% tap amplitudes in this case.
randomRatioForAmplitude = 0;
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