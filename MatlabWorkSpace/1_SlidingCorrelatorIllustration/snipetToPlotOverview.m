% SNIPETTOPLOTOVERVIEW A helper script to generate the overview figure.
%
% Yaguang Zhang, Purdue, 07/17/2021

% Pattern start points.
sigTxStartTs = 1:(simTimeLength*NUM_OF_TX_SIG_CYCLE_TO_SHOW);
sigTxStartTs = sigTxStartTs(mod(sigTxStartTs, sigTxCycleLen)==1);
sigRxStartTs = 1:(simTimeLength*NUM_OF_TX_SIG_CYCLE_TO_SHOW);
sigRxStartTs = sigRxStartTs(mod(sigRxStartTs, sigRxCycleLen)==1);

sigTxToShow = repmat(sigTx, [NUM_OF_TX_SIG_CYCLE_TO_SHOW, 1]);
sigRxToShow = repmat(sigRx, [NUM_OF_TX_SIG_CYCLE_TO_SHOW, 1]);
sliCorToShow = repmat(sliCor, [NUM_OF_TX_SIG_CYCLE_TO_SHOW, 1]);

curFig = figure;
hSPTx = subplot(3,1,1); hold on; grid on; grid minor;
plot(sigTxToShow, 'Color', 'g', 'LineWidth', normLineW);
xticks([]); axis tight; ylabel('TX signal');
curAxis = axis;
plot([sigTxStartTs; sigTxStartTs], ...
    repmat(curAxis(3:4)', [1, length(sigTxStartTs)]), ...
    startLinePattern, 'Color', startLineC, 'LineWidth', thickLineW);
hSPTx.YLimMode = 'manual';
ax = gca; 
ax.FontSize = subPlotFontSize;
ax.FontWeight = 'bold';

hSPRx = subplot(3,1,2); hold on; grid on; grid minor;
plot(sigRxToShow, 'Color', 'b', 'LineWidth', normLineW);
xticks([]); axis tight; ylabel('RX signal');
curAxis = axis;
plot([sigRxStartTs; sigRxStartTs], ...
    repmat(curAxis(3:4)', [1, length(sigRxStartTs)]), ...
    startLinePattern, 'Color', startLineC, 'LineWidth', thickLineW);
hSPRx.YLimMode = 'manual';
ax = gca; 
ax.FontSize = subPlotFontSize;
ax.FontWeight = 'bold';

hSPCo = subplot(3,1,3); hold on; grid on; grid minor;
plot(sliCorToShow, 'Color', 'k', 'LineWidth', normLineW);
xticks([]); axis tight; ylabel({'Windowed'; 'Correlation'});
xlabel('Time'); hold on; grid on; grid minor;
hSPCo.YLimMode = 'manual';
ax = gca; 
ax.FontSize = subPlotFontSize;
ax.FontWeight = 'bold';

sgtitle(['N = ', num2str(N)]);

linkaxes([hSPTx hSPRx hSPCo], 'x');
% EOF