function do_EEGLABPlot_masked(ersp, erspboot, times, freqs, ersplim, MarkTimes, p_values, titleString, plotThresholded, varargin)

% Plot the ERSP (Event-Related Spectral Perturbation)
if plotThresholded == 1
    % Zero out periods that do not reach significance
    for ii = 1:length(erspboot)
        thresholdedersp(ii,:) = ersp(ii,:).*(ersp(ii,:)<(erspboot(ii,1)) | ersp(ii,:)>(erspboot(ii,2)));
    end
    imagesc(times, freqs, thresholdedersp(:,:), ersplim);
else
    imagesc(times, freqs, ersp(:,:), ersplim);
end

set(gca, 'YDir', 'normal', 'YScale', 'log');
xlabel('Time (ms)');
ylabel('Frequency (Hz)');

% Set frequency ticks
tick = [4, 8, 15, 25, 50, 100, 150]; % Customize as needed
set(gca, 'YTick', tick);
hold on;

% Plot mark times as dotted lines
for mt = MarkTimes(:)'
    plot([mt mt], [0.01, freqs(end)], '--k', 'LineWidth', 2);
end

% Set title
title(titleString);

% Set y-limits to min/max tick
if nargin >= 10
    ylim(varargin{1});
else
    ylim([tick(1), tick(end)]);
end

% Add colorbar
colorbar();

% Create and overlay the semi-transparent mask based on p-values
mask = p_values > 0.05; % True where p-values are greater than 0.05
h = imagesc(times, freqs, double(mask)); % Convert logical to double for imagesc
set(h, 'AlphaData', mask * 1); % Apply semi-transparency to non-significant areas
set(gca, 'YDir', 'normal'); % Ensure y-direction is normal after overlaying the mask

end
