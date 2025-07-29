function plotPSDGroupComp(TM, NPX, opt)
% plotPSDGroupComp Compare real PSD between vigilance states across two EpyData instances.
%
% Inputs:
%   - TM:    EpyData object representing Telemetry (order matters)
%   - NPX:   EpyData object representing Neuropixels (order matters)
%   - opt:   Optional structure with fields:
%              • ratTM:   IDs of Telemetry rats to include for concatenated PSD
%              • ratNPX:  IDs of Neuropixels rats to include for concatenated PSD
%              • binCount: number of time bins for PSD estimation
%              • range:   analysis time range (default: 0 = full session)
%
% Output:
%   First Figure
%   - Generates a figure with 3 subplots:
%       1. Sleep vs Wake PSD for Telemetry rats
%       2. Sleep vs Wake PSD for Neuropixels rats
%       3. Sleep vs Wake PSD for the pooled population (Telemetry + NPX)
%
%   - Each subplot shows log-period vs PSD curves (semplotPietro), with comparisons
%     across vigilance states and spectral signatures.
%
%   - PSDs are computed per rat, averaged per group and condition,
%     then visualized in log-scale over period (1/f in days).
%
%   Seconde Figure
%   - Displays statistical comparison of PSD values over the all population between vigilance states
%     (Sleep vs Wake) across time-frequency bins using unpaired Wilcoxon tests.
%   - Shows p-value maps (bin-by-bin) for each population (TM and NPX),
%     highlighting significant differences (e.g., p < 0.05).n.
%
%   Third Figure
%   - Same as the second one but ONLY with the TM population
%
% Dependencies :
%   - myBoxPlox
%   - computePSD


arguments
    TM (1,1) EpyData
    NPX (1,1) EpyData
    opt.ratTM(1,:) double = TM.ratID                      % List of Telemetry rat IDs
    opt.ratNPX(1,:) double = NPX.ratID                    % List of Neuropixel rat IDs
    opt.binCount (1,1) double = 3600                      % Number of time bins for PSD
end

states = ["sleep","wake"];                                % Define vigilance states
nStates = numel(states);
groups = ["Telemetry", "Neuropixel", "All"];              % Comparison groups
nGroups = numel(groups);

% --- Retrieve PSDs for each group (TM, NPX) and each state ---
for g = 1:nGroups-1                                       % Skip "All" for now
    group = groups(g);
    for s = 1:nStates
        state = states(s);
        switch group
            case "Telemetry"
                [~, PSDReal.(group).(state)] = TM.computePSD("rat", opt.ratTM, ...
                    "state", state, "range",TM.rangeCutIS, "binCount", opt.binCount);
            case "Neuropixel"
                [~, PSDReal.(group).(state)] = NPX.computePSD("rat", opt.ratNPX, ...
                    "state", state, "range", NPX.rangeCutIS, "binCount", opt.binCount);
        end
        % Concatenate real PSDs across rats into matrix: freq x rats
        pxxReal.(group).(state) = cell2mat(cellfun(@(x) x.pxx, ... % ./max(x.pxx), ...% each PSD not normalized
            PSDReal.(group).(state), 'UniformOutput', false)');
    end

    % Find the max between pxx sleep and pxx wake for each rats to find the
    % normalisation factors
    pxxReal.(group).normFactors = max(max(pxxReal.(group).sleep), max(pxxReal.(group).wake));

    % Extract frequency vector from one PSD (they're all the same for a given group)
    pxxReal.(group).f = PSDReal.(group).(state){1}.f;

end

% The f grid to compute the PSD of TM is a fortiori strictly included in
% the f grid of NPX. We then limit the comparison to the smallest fGrid. To
% understand this comparison limitation requirement check the method
% computePSD and the construction of the FGrid for PSD
sharedFreq = pxxReal.Neuropixel.f; 

% --- Pool all rats across both modalities into "All" group ---
for s = 1:nStates
    state = states(s);
    pxxTM = pxxReal.Telemetry.(state);
    pxxNPX = pxxReal.Neuropixel.(state);
    pxxReal.All.(state) = [pxxTM(end - numel(sharedFreq)+1:end,:),... % Share of the "highest" frequency     
        pxxNPX]; 
    pxxReal.All.f = sharedFreq;
    pxxReal.All.normFactors = [pxxReal.Telemetry.normFactors, pxxReal.Neuropixel.normFactors];
end

% --- Plot each group’s PSD: sleep vs wake ---
figure('Units', 'pixels', 'Position', [0, 0, 2000, 1200]);
for g = 1:nGroups
    subplot(3,1,g);

    % Plot sleep PSDs (log-period scale)
    semplotPietro(log(1./(pxxReal.(groups(g)).f*86400)), pxxReal.(groups(g)).sleep./pxxReal.(groups(g)).normFactors, ...
        [1.00, 0.00, 0.60], 'mode', 'std', 'legend', 'Sleep PSD');

    hold on;

    % Plot wake PSDs
    semplotPietro(log(1./(pxxReal.(groups(g)).f*86400)), pxxReal.(groups(g)).wake./pxxReal.(groups(g)).normFactors, ...
        [0.20, 0.40, 1.00], 'mode', 'std', 'legend', 'Wake PSD');

    % Axis labels and title
    xlabel('Time (days, log scale)');
    ylabel('Normalized PSD (a.u.) ');
    title(sprintf('Normalized PSD for group %s', groups(g)));

    % Legend and aesthetics
    legend('Location', 'northwest');
    grid on;
    ylim([0 1]);
    xticks([log(1/2), log(1), log(3), log(7), log(10), log(15),log(21)]);
    xticklabels({'12h', '24h', '3d', '7d', '10d', '15d', '21d'});
    switch groups(g)
        case "Telemetry"
            xlim([-1, log(21)]); 
        otherwise % NPX and ALL comparison
            xlim([-1, log(15)]);
    end


    % --- Compute Wilcoxon rank-sum (unpaired) for each frequency bin ---
    % Frequency limit to do stat
    freqBounds = [pxxReal.(groups(g)).f(1), 1/(3*24*3600)]; % [min of the f grid (f15 or f21 days), f3 days]
    frestrict = Restrict(sharedFreq,freqBounds);
    ps = nan(size(frestrict));                                      % p-values vector
    for i = 1:length(frestrict)
        ps(i) = ranksum(pxxReal.(groups(g)).sleep(i, :), ...
            pxxReal.(groups(g)).wake(i, :));    % Unpaired test
    end

    % Highlight significant bins (e.g., p < 0.05) ---
    hold on;
    sigBins = ps < 0.05;
    yMax = ylim;
    scatter(log(1./(frestrict(sigBins)*86400)), ...
        yMax(2)*0.98*ones(1,sum(sigBins)), ...
        10, 'k', 'filled', DisplayName="Significative frequencies");  % Add small black dots at top
end

sgtitle('PSD by Groups');  % Global title


% --- Save figure ---
figName = "Sleep_VS_WAKE_PSD_byGroup";
savePath = fullfile("/mnt/hubel-data-103/Corentin/Matlab/Figures_PSD", figName);
exportgraphics(gcf, savePath + ".png", 'Resolution', 300);
    % saveas(gcf, savePath + ".png");
    % saveas(gcf, savePath + ".tif");

%% Comparaison

dataTM = cell2mat(cellfun(@(x) [x.incrFromShufSleep(1), ... % Peak increase of first most important peak - sleep
                     x.incrFromShufWake(1), ... % Peak increase of first most important peak - wake
                     sum(x.incrFromShufSleep), ... % Cumulative peak increase of significative peaks - sleep
                     sum(x.incrFromShufWake), ... % Cumulative peak increase of significative peaks - wake
                     x.areaFromShufSleep(1), ... % Area increase of first most important peak - sleep
                     x.areaFromShufWake(1), ... % Area increase of first most important peak - wake
                     sum(x.areaFromShufSleep), ... % Cumulative area increase of significative peaks - sleep
                     sum(x.areaFromShufWake)], ... % Cumulative area increase of significative peaks - wake
                     TM.multidienCycl(cellfun(@(x) ~isempty(x) && ~isempty(x.f0),TM.multidienCycl)), 'UniformOutput',false));

dataNPX = cell2mat(cellfun(@(x) [x.incrFromShufSleep(1), ... % Peak increase of first most important peak - sleep
                     x.incrFromShufWake(1), ... % Peak increase of first most important peak - wake
                     sum(x.incrFromShufSleep), ... % Cumulative peak increase of significative peaks - sleep
                     sum(x.incrFromShufWake), ... % Cumulative peak increase of significative peaks - wake
                     x.areaFromShufSleep(1), ... % Area increase of first most important peak - sleep
                     x.areaFromShufWake(1), ... % Area increase of first most important peak - wake
                     sum(x.areaFromShufSleep), ... % Cumulative area increase of significative peaks - sleep
                     sum(x.areaFromShufWake)], ... % Cumulative area increase of significative peaks - wake
                     NPX.multidienCycl(cellfun(@(x) ~isempty(x) && ~isempty(x.f0),NPX.multidienCycl)), 'UniformOutput',false));

xlabels = {"sleep", "wake"};
titles = {"Main peak increase from surrogates",...
        "Cumulative peak increase from surrogates",...
        "Main area increase from surrogates",...
        "Cumulative area increase from surrogates"};
ylabels = {"Distance (a.u.)",...
        "Cumulative distance (a.u.)",...
        "Area (a.u.)",...
        "Cumulative area (a.u.)"};

%% Both population
figure('Name', 'TM & NPX - Increase normalized power form shuffle sleep VS wake', ...
        'Units', 'pixels', 'Position', [0, 0, 2400, 1600], ...
        'Visible','on');
hold on;


for d = 1:floor(size(dataTM,2)/2)
    subplot(floor(size(dataTM,2)/4), floor(size(dataTM,2)/4),d)
    metricSleep = [dataTM(:,1+2*(d-1)); dataNPX(:,1+2*(d-1))];
    metricWake = [dataTM(:,2+2*(d-1)); dataNPX(:,2+2*(d-1))];
    
    % BoxPlot
    group = [ones(size(metricSleep)); 2*ones(size(metricWake))];
    myBoxPlot([metricSleep; metricWake], group, ...
        'xtlabels', xlabels, ...
        'box', 2, 'boxcolors', 'w', 'scatter', 0, ...
        'outliers', 0, ...
        'jitter', 0, 'violin', 'none', 'colors', [0 0 0]);

    % Significance of boxplot
    [p,h]= signrank(metricSleep,metricWake);

    % Coordinates for significance bar
    yMax = max([metricSleep; metricWake]) * 1.1;
    lineX = [1, 2];
    lineY = [yMax, yMax];
    textY = yMax * 1.05;

    hold on;
    plot(lineX, lineY, 'k-', 'LineWidth', 1.5, 'HandleVisibility','off');

    if h == 0 || isnan(p)
        text(mean(lineX), textY, 'n.s.', ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k');
    else
        % Compute number of stars based on p-value
        nStars = floor(-log10(p));
        nStars = max(nStars, 1);        % At least one star
        nStars = min(nStars, 4);        % Avoid too many stars
        stars = repmat('*', 1, nStars);

        % Display stars and p-value
        text(mean(lineX), textY, sprintf('%s (p=%.1e)', stars, p), ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k');
    end

    % Scatter of data
    groupTM = [ones(size(dataTM(:,1+2*(d-1)))); 2*ones(size(dataTM(:,2+2*(d-1))))];
    scatterTM = myBoxPlot([dataTM(:,1+2*(d-1));dataTM(:,2+2*(d-1))], groupTM, ...
        'xtlabels', xlabels, ...
        'box', 0, 'scatter', 2, ...
        'jitter', 0, 'withinlines', 1, 'violin', 'none', 'colors',[0 0 1]);
   
    groupNPX = [ones(size(dataNPX(:,1+2*(d-1)))); 2*ones(size(dataNPX(:,2+2*(d-1))))];
    scatterNPX = myBoxPlot([dataNPX(:,1+2*(d-1));dataNPX(:,2+2*(d-1))], groupNPX, ...
        'xtlabels', xlabels, ...
        'box', 0, 'scatter', 2, ...
        'jitter', 0, 'withinlines', 1, 'violin', 'none', 'colors',[1 0 0]);

    legend([scatterTM.sc(1),scatterNPX.sc(1)], ["TM", "NPX"])
    sgtitle('TM & NPX - Increase normalized power form shuffle sleep VS wake')
    title(titles(d))
    ylabel(ylabels(d))
    ylim([min([metricSleep; metricWake]), yMax*1.1])
end
%% TM population
figure('Name', 'TM - Increase normalized power form shuffle sleep VS wake', ...
        'Units', 'pixels', 'Position', [0, 0, 2400, 1600], ...
        'Visible','on');
hold on;

for d = 1:floor(size(dataTM,2)/2)
    subplot(floor(size(dataTM,2)/4), floor(size(dataTM,2)/4),d)
    metricSleep = [dataTM(:,1+2*(d-1))];
    metricWake = [dataTM(:,2+2*(d-1))];
    
    % BoxPlot
    group = [ones(size(metricSleep)); 2*ones(size(metricWake))];
    myBoxPlot([metricSleep; metricWake], group, ...
        'xtlabels', xlabels, ...
        'box', 2, 'boxcolors', 'w', 'scatter', 0, ...
        'jitter', 0, 'violin', 'none', 'colors', [0 0 0]);

    % Significance of boxplot
    [p, h]= signrank(metricSleep,metricWake);

    % Coordinates for significance bar
    yMax = max([metricSleep; metricWake]) * 1.1;
    lineX = [1, 2];
    lineY = [yMax, yMax];
    textY = yMax * 1.05;

    hold on;
    plot(lineX, lineY, 'k-', 'LineWidth', 1.5, 'HandleVisibility','off');

    if h == 0 || isnan(p)
        text(mean(lineX), textY, 'n.s.', ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k');
    else
        % Compute number of stars based on p-value
        nStars = floor(-log10(p));
        nStars = max(nStars, 1);        % At least one star
        nStars = min(nStars, 4);        % Avoid too many stars
        stars = repmat('*', 1, nStars);

        % Display stars and p-value
        text(mean(lineX), textY, sprintf('%s (p=%.1e)', stars, p), ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'k');
    end

    % Scatter of data
    groupTM = [ones(size(dataTM(:,1+2*(d-1)))); 2*ones(size(dataTM(:,2+2*(d-1))))];
    scatterTM = myBoxPlot([dataTM(:,1+2*(d-1));dataTM(:,2+2*(d-1))], groupTM, ...
        'xtlabels', xlabels, ...
        'box', 0, 'scatter', 2, ...
        'jitter', 0, 'withinlines', 1, 'violin', 'none', 'colors',[0 0 1]);
  
    legend([scatterTM.sc(1)], "TM")
    title(titles(d))
    sgtitle('TM - Increase normalized power form shuffle sleep VS wake')
    ylabel(ylabels(d))
    ylim([min([metricSleep; metricWake]), yMax*1.1])
end

end