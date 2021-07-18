% SNIPETTOPLOTOVERVIEW A helper script to generate multiple illustrations.
%
% Yaguang Zhang, Purdue, 07/17/2021

for N = 2.^(2:8)-1
    FLAG_GEN_ANIME = N<128;
    illustrateSlidingCorrelator;
end

% EOF