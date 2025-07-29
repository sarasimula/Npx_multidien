function plotPSDPerRat(this, opt)
% plotPSDPerRat - Plots IS time series and PSDs (real + surrogates) per vigilance state and per rat.
%                 Generates one figure per rat, combining time series and PSD visualization.
%
% SYNTAX:
%   this = plotPSDPerRat(this)
%   this = plotPSDPerRat(this, opt)
%
% INPUTS:
%   this        : (1,1) EpyData object, containing fields:
%                   - PSDreal / PSDsurro : cell of structs with fields `.f` and `.pxx` per state
%                   - multidienCycl      : cell with fields `freqISFiltered`, `f0`, `fWidth`, `sigBandEdges`
%                   - seizuresTimestamps : cell with seizure event data per rat
%                   - refDatesTimestamps : datetime ranges used for temporal alignment
%                   - ratID              : vector of rat IDs
%                   - samplingFreq       : double, modality-specific sampling frequency
%   opt         : struct with optional parameters:
%                   - rat      : vector of rat IDs to process (default: all rats in `this.ratID`)
%                   - binCount : bin size (in seconds) for binning IS counts (default: 3600)
%
% OUTPUT:
%   this        : Updated EpyData object (no new fields added, but plotting may rely on precomputed values)
%
% SIDE EFFECTS:
%   - Creates and saves figures (PNG format) under hardcoded path `/mnt/hubel-data-103/Corentin/Matlab/Figures_PSD`
%   - Figures are saved using rat ID and inferred modality (NPX or TM)
%
% DEPENDENCIES (called functions/methods):
%   - this.getFrequency        % method of EpyData
%   - signalProcessing         % smooths and interpolates IS time series using harmonics
%   - overlaySeizures          % overlays seizure information on time series plot
%   - semplotPietro                 % custom plotting of surrogates (mean ± SEM or shaded)
%
% NOTES:
%   - PSDs are normalized by their individual max before plotting
%   - x-axis in PSD plots uses log(period in days)
%   - Highlighting includes statistically significant bands
%
% AUTHOR:
%   Corentin L., 2025
%
% SEE ALSO:
%   computePSD, analyzePSD, getFrequency

arguments
    this (1,1) EpyData
    opt.rat(1,:) double = this.ratID
    opt.binCount (1,1) double = 3600
end

states = ["all", "sleep", "wake"];
nStates = numel(states);
fMax = 1 / opt.binCount;

for i = 1:numel(opt.rat)
    ratID = opt.rat(i);
    idx = find(this.ratID == ratID);
    if isempty(idx), continue; end
    
    figure('Name', sprintf('Rat %d - Temporal Evolution and PSD by State', ratID), ...
        'Units', 'pixels', 'Position', [0, 0, 2400, 1600], ...
        'Visible','off');

    axLeft = gobjects(nStates,1);  % Handle the y axis for IS time series (figures de gauche)
    axRight = gobjects(nStates,1); % Handle the y axis for PSD (figures de droite)
    [YminPSD, YmaxPSD, YminCount, YmaxCount] = deal(-inf);

    for s = 1:nStates
        state = states(s);

        % fReal = PSDReal{1}.f;
        fReal = this.PSDreal{idx}.(state).f;
        fSurro = this.PSDsurro{idx}.(state).f;
        idxReal = fReal <= fMax;
        idxSurro = fSurro <= fMax;

        % pxxReal = PSDReal{1}.pxx/max(PSDReal{1}.pxx); 
        pxxReal = this.PSDreal{idx}.(state).pxx(idxReal)./max(this.PSDreal{idx}.(state).pxx(idxReal)); % Normalized real PSD by its max 
        pxxSurro = cell2mat(cellfun(@(x) x(idxSurro)./max(x(idxSurro)), this.PSDsurro{idx}.(state).pxx, 'UniformOutput', false)'); % Matrix of normalized surrogates PSD by its respective max 

        % Keep track for ylim axis
        YmaxPSD = max([YmaxPSD, max(pxxReal(idxReal))]);
        YminPSD = min([YminPSD, min(pxxReal(idxReal))]);

        % -- Plot: IS frequency over time --
        freq = this.getFrequency("rat", ratID, "nature", "IS", ...
            "state", state, "binCount", opt.binCount).freq{1}'; % * opt.binCount;
        t = (0:numel(freq)-1) / 24; % time in days

        freqProcessed = signalProcessing(freq, nHarmonics=20);

        % Keep track for ylim axis
        YmaxCount = max([YmaxCount, max(freqProcessed)]);
        YminCount = min([YminCount, min(freqProcessed)]);

        axLeft(s) = subplot(nStates, 2, 2*s - 1);
        hold on;
        
        plot(t, freq, 'Color', [0 0 1 0.3], 'LineWidth', 1.5, 'DisplayName', 'Raw - used for PSD');
        plot(t, freqProcessed, 'Color', [1 0.5 0], 'LineWidth', 2.5, 'DisplayName', 'Interpoled and Smoothed - not used for PSD');
        
        % Sinus, seizure and state overlay (only in 'all' state)
        if s == 1 
            % Seizures
            if ~isempty(this.seizuresTimestamps{idx})
                seiz = this.seizuresTimestamps{idx};
                tEnd = seconds(this.refDatesTimestamps(idx,2) - this.refDatesTimestamps(idx,1));

                overlaySeizures(seiz, tEnd, numberHours = 6, ...
                    eventDisp = true, cumuCountDisp = true)
            end

            % State
            % overlayStates(this,idx)
        end


        % Multidien signal new method
        if ~isempty(this.multidienCycl{idx,1}.freqISFiltered)
            plot(t, this.multidienCycl{idx,1}.freqISFiltered.(state), 'Color', [0.1 0.2 0.1], 'LineWidth', 2.5, 'DisplayName', 'Multidien Rythm');
        end
            % % Old Method : Sinus fitted when there was only one
            % omega = 2 * pi * this.multidienCycl{idx,1}.f0 * (24*3600); % in d⁻¹
            % A = this.multidienCycl{idx,1}.A;        % amplitude
            % phi_t0 = this.multidienCycl{idx,1}.phi_t0;     % phase
            % fit_sin = A * cos(omega*t + phi_t0) + mean(freqFiltered);
            % plot(t, fit_sin, 'Color', [0.1 0.2 0.1], 'LineWidth', 2.5, 'DisplayName', 'Best LF sinusoïd')
 
        xlabel('Time (days)');
        xlim([0,seconds(this.refDatesTimestamps(idx,2)- this.refDatesTimestamps(idx,1))/(24*3600)])
        ylabel(sprintf('Normalized IS count per bin of %d s', opt.binCount));
        title(sprintf('IS Time Series - %s', state));
        legend;
        grid on;

        % -- Plot: PSD spectrum with surrogates and significant bands --
        axRight(s) = subplot(nStates, 2, 2*s);
        plot(log(1./(fReal(idxReal)*86400)),... % Log scale for x axis
            pxxReal(idxReal), ... % Already normalized in the beggining
            'LineWidth', 2, 'DisplayName','Real', 'Color', [1.00, 0.80, 0.00]);
        hold on;
        semplotPietro(log(1./(fSurro(idxSurro)*86400)), ... % Log scale for x axis
            pxxSurro, ...  % Already normalized in the beggining
            [1.00, 0.30, 0.30], 'mode', 'sem', 'legend','Surrogate');

        fMultidien = this.multidienCycl{idx,1}.f0;
        fWidth = this.multidienCycl{idx,1}.fWidth;
        if ~isempty(fMultidien)
            for l = 1:numel(fMultidien)
                % Plot in green the multidien frequency
                xline(log(1/(fMultidien(l)*86400)), 'g--', 'LineWidth', 0.75, 'HandleVisibility', 'off'); 
                % Plot in black the frequencies used to define the width of the band pass of the filter
                xline(log(1/(fWidth(l,:)*86400)), 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off');
                hold on;
            end
        end
        
        xlabel('Period (days, log scale)');
        ylabel('PSD (a.u.)');
        title(sprintf('PSD - %s', state));
        legend('Location','northwest');
        grid on;

        switch this.samplingFreq
            case 512 % TM
                xlim([-1, log(21)]); % until 21 days
            case 30000 % NPX
                xlim([-1, log(15)]); % until 15 days
        end
        
        xticks([log(1/2), log(1), log(3), log(7), log(10), log(15), log(21)]);
        xticklabels({'12h', '24h', '3d', '7d', '10d', '15d', '21d'});

        % Highlight statistically significant bands (above surrogate)
        sigBandEdges = this.multidienCycl{idx,1}.sigBandEdges;
        if ~isempty(sigBandEdges)
            yMax = max(pxxReal(idxReal)./max(pxxReal(idxReal))) * 1.02;
            for b = 1:size(sigBandEdges,1)
                fL = sigBandEdges(b,1);
                fH = sigBandEdges(b,2);
                xL = log(1 / (fH * 86400));
                xH = log(1 / (fL * 86400));
                plot([xL, xH], [yMax, yMax], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
            end
        end
    end

    % Compute adjusted edges
    ylimCount = [YminCount - 0.25*abs(YminCount), YmaxCount + 0.05*abs(YmaxCount)];
    ylimPSD   = [YminPSD,   YmaxPSD   + 0.08*abs(YmaxPSD)];

    % Apply limits
    for s = 1:nStates
        ylim(axLeft(s), ylimCount);
        ylim(axRight(s), ylimPSD);
        xlim(axLeft(s), [0,this.rangeCutIS(idx)/(3600*24)])
    end

    sgtitle(sprintf('Rat %d – Temporal Evolution and PSD by State', ratID));
    set(gcf, 'DefaultTextInterpreter', 'none');

    % -- Modality identification based on sampling frequency --
    if this.samplingFreq == 30000
        modality = "NPX"; % Neuropixels
    elseif this.samplingFreq == 512
        modality = "TM";  % Telemetry
    end

    % -- Save the figure with informative filename --
    figName = sprintf("v10-%02d_%s_PSDandTimeSeries-raw", ratID, modality);
    savePath = fullfile("/mnt/hubel-data-103/Corentin/Matlab/Figures_PSD", figName);
    exportgraphics(gcf, savePath + ".png", 'Resolution', 300);
end
end


%%        %        %         %        %    HELP FUNCTIONS    %         %        %        %         % 

function overlayStates(this, idx)
% overlayStates  Overlay sleep and wake periods on current plot.
%
%   overlayStates(this, idx) shades the regions corresponding to sleep and
%   wake periods from the object's sleep and wake properties. Sleep and wake intervals are displayed as semi-
%   transparent patches on the current axes (gca).
%
%   Inputs:
%       this : Object containing `sleep` and `wake` properties,
%              each a cell array of Nx2 matrices of time intervals.
%       idx  : Index selecting which session's states to overlay.
%
%   Notes:
%       - Sleep and wake periods are defined by their start and end times.
%       - This function modifies the current axes by adding shaded regions.
%
%   See also: patch, gca

sleep = this.sleepInt{idx};
wake = this.wakeInt{idx};
yLimits = get(gca, 'YLim');
for i = 1:numel(sleep(:,1))-1 % Avoid log(0)
     x0 = sleep(i+1,1)/(3600*24); % In days
     x1 = sleep(i+1,2)/(3600*24); % In days
     vis = 'on'; lbl = 'sleep intervals';
     if i > 1
         vis = 'off'; lbl = '';
     end
     patch(gca, [x0 x1 x1 x0], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
                [0.7,0.7,0] , 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', vis, 'DisplayName', lbl);
end

for j = 1:numel(wake(:,1))-1 % Avoid log(0)
     x0 = wake(j+1,1)/(3600*24); % In days
     x1 = wake(j+1,2)/(3600*24); % In days
     vis = 'on'; lbl = 'wake intervals';
     if j > 1
         vis = 'off'; lbl = '';
     end
     patch(gca, [x0 x1 x1 x0], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
                [0,1,0] , 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', vis, 'DisplayName', lbl);
end

end