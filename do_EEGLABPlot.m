% use this for plotting EEGLAB wavelet data, as the way it does it isn't
% very nice (e.g. can't plot log scale whilst including markers)

function do_EEGLABPlot(ersp,erspboot,times,freqs,ersplim,MarkTimes,titleString,plotThresholded,varargin)

% Plot the figure (whether thresholded or not)
% colormap(jet(256));
if plotThresholded == 1
    % THIS ZEROES OUT ANY PERIODS THAT DO NOT REACH
    % SIGNIFICANCE (AS DETERMINED BY ERSPBOOT)
    for ii = 1:length(erspboot)
        thresholdedersp(ii,:) = ersp(ii,:).*(ersp(ii,:)<(erspboot(ii,1))|ersp(ii,:)>(erspboot(ii,2)));
    end
    imagesc(times,freqs,thresholdedersp(:,:), ersplim)
else
    imagesc(times,freqs,ersp(:,:), ersplim)
end

set(gca,'YDir','normal','YScale', 'log')
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
tick = [2 4 8 15 25 50 100 150];
tick = [4 8 15 25 50 100 150];
set(gca,'YTick',tick);  
hold on

% plot mark times dotted lines (from y=0.01 because of log scale)
for mt = MarkTimes(:)'
    plot([mt mt],[0.01 freqs(end)],'--k','LineWidth',2);
end

% set main title
title(titleString)

% set y-lim to min/max tick
if nargin == 9
    ylim(varargin{1})
else
    ylim([tick(1) tick(end)])
end
% add colourbar
h(2) = gca;
h3 = colorbar();
% h(3) = cbar('vert'); % ERSP colorbar axes
% title('ERSP(dB)')