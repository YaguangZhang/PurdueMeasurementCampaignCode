function [ h ] = plotPwmPulse( PWM )
%PLOTPWMPULSE Plot a sequence in the form of a pulse width modulation
%signal.
%  
% Ref: 
%   https://www.mathworks.com/matlabcentral/answers/374185-plot-rectangular-square-wave
%
% Yaguang Zhang, Purdue, 03/14/2019

PWM = PWM(:)';
t_V = 1:numel(PWM);
idx = diff(PWM)~=0;
X = repmat(t_V([true,idx]),2,1);
Y = repmat(PWM([true,idx]),2,1);
X = X(2:end);
X = [X X(end)+1];
Y = Y(1:end-1);
Y = [Y Y(end)];
h = plot(X,Y);

axis tight;
curAxis = axis;
deltaYToExtendEachSide = (curAxis(4)-curAxis(3))/10;
axis([curAxis(1:2) ...
    curAxis(3)-deltaYToExtendEachSide curAxis(4)+deltaYToExtendEachSide]);
grid minor;

end
% EOF