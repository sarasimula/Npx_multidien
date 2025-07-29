function overlaySeizures(seiz, tEnd, opt)
% overlaySeizures - Overlays seizure events and seizure count on an existing time series plot.
%
% This function adds vertical lines at the midpoints of seizure events and optionally
% overlays a histogram of seizure counts over time using a secondary Y-axis.
%
% Inputs:
%   seiz            : (n,2) double, start and end timestamps of seizures [in seconds]
%   tEnd            : (1,1) double, total session duration [in seconds]
%   opt.numberHours : (1,1) double, time bin size for seizure count display [in hours], default = 6
%   opt.eventDisp   : (1,1) logical, whether to display vertical lines for seizure events, default = true
%   opt.cumuCountDisp : (1,1) logical, whether to display cumulative seizure count, default = true

    arguments
        seiz (:,2) double
        tEnd (1,1) double {mustBePositive}
        opt.numberHours (1,1) double {mustBePositive} = 6
        opt.eventDisp (1,1) logical = true
        opt.cumuCountDisp (1,1) logical = true
    end

    % Compute seizure midpoints (in days)
    seizMid = mean(seiz, 2) / (3600 * 24);

    % Optionally overlay vertical dashed lines at each seizure midpoint
    if opt.eventDisp
        % First line different for the legend
        xline(seizMid(1), '--r', 'LineWidth', 0.5, 'Alpha', 0.5, 'DisplayName', 'Seizures');
        xline(seizMid(2:end), '--r', 'LineWidth', 0.5, 'Alpha', 0.5, 'HandleVisibility','off');
    end
   
    % Optionally overlay cumulative seizure count using right Y-axis
    if opt.cumuCountDisp
        binSeizures = opt.numberHours * 3600; % Bin width in seconds
        edges = 0:binSeizures:tEnd;
        countSeiz = histcounts(seizMid, edges);
        binCenters = (edges(1:end-1) + edges(2:end)) / 2 / (3600*24); % Convert to days

        yyaxis right;
        plot(binCenters, countSeiz, '-s', 'Color', [1 0 0], ...
            'LineStyle', '-', 'LineWidth', 1.5, 'HandleVisibility','off');
        ylabel('Seizure count');
        
        ax = gca;
        ax.YColor = [1 0 0];
        ylim([0 30]); % Optional: adjust to fit your data
        yyaxis left;
    end
end