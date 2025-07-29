function plotGlobalEventFrequencyDistribution(TM, NPX, ratTM, ratNPX, opt)
% plotGlobalEventFrequencyDistribution - This function compares the distribution of interictal spikes (IS) and seizure
% frequencies between two experimental groups: one recorded using telemetry (TM)
% and the other using Neuropixels (NPX). It generates a multi-panel violin plot
% summarizing event rates across different vigilance states (all, sleep, wake)
% and saves the figure in both .svg and .tif formats.
%
% SYNTAX : 
%   plotGlobalEventFrequencyDistribution(TM, NPX, ratTM, ratNPX, saveDir)
%
% INPUTS:
%   TM      : An instance of the EpyData class containing telemetry data.
%             Must support the method getFrequency(...).
%
%   NPX     : An instance of the EpyData class containing Neuropixels data.
%             Must support the method getFrequency(...).
%
%   ratTM   : cell of telemetry rat IDs to include in the comparison.
%             {ratNoEmptySeizTM, ratNoEmptyStateTM}
%
%   ratNPX  : cell of Neuropixels rat IDs to include in the comparison.
%             {ratNoEmptySeizNPX, ratNoEmptyStateNPX}
%
% FUNCTIONALITY:
%   - Generates violin plots of average seizure frequency (1/day) and IS frequency (1/h)
%     across rats for both groups.
%   - Subdivides IS analysis by vigilance state (all, sleep, wake).
%   - Adds a final subplot comparing sleep vs wake IS rates across all rats pooled.
%   - Applies standard plot formatting and uses external functions:
%       * plotComparativeViolin(dataCell, xLabels)
%       * beautifyPlot()
%   - Saves the figure under the name "GlobalDistributionOverPopulations".
%
% DEPENDENCIES:
%   - The EpyData objects must implement:
%       * getFrequency("rat", ratList, "nature", natureStr, "state", stateStr, "binCount", binVal)
%           which returns a struct with an 'avg' field (average event rate per rat).
%
% EXAMPLE USAGE:
%   plotGlobalEventFrequencyDistribution(TM, NPX, ratNoEmptySeizTM, ratNoEmptySeizNPX, ...
%       '/mnt/hubel-data-103/Corentin/Matlab/Figures_PSD');
%
% AUTHOR:
%   Corentin - Data analysis for multiscale epilepsy biomarkers.
%
% VERSION:
%   July 2025


arguments
    TM (1,1) EpyData
    NPX (1,1) EpyData
    ratTM (1,:) cell
    ratNPX (1,:) cell
    opt.binCountIS (1,1) double = 3600
    opt.binCountSeiz (1,1) double = 86400
end 

% Remove rat 800 from NPX
ratNPX{1} = setdiff(ratNPX{1}, 800);
ratNPX{2} = setdiff(ratNPX{2}, 800);

% Parameters
nature = ["IS","seizures"];
state = ["all","sleep","wake"];
xlabels = {'Telemetry','Neuropixel'};

% Figure setup
height_cm = 100;
width_cm = 1.5 * height_cm;
fig = figure('Units', 'centimeters', 'Position', [0 0 width_cm height_cm]);
nPlots = (numel(nature) - 1) + numel(state); 
nCols = 3;
nRows = ceil(nPlots / nCols);
plotIdx = 1;

% Loop over event types
for i = 1:numel(nature)
    if nature(i) == "seizures"
        subplot(nRows, nCols, plotIdx);
        d1 = TM.getFrequency("rat", ratTM{1}, "nature", "seizures", "state", "all", "binCount", opt.binCountSeiz).avg;
        d2 = NPX.getFrequency("rat", ratNPX{1} , "nature", "seizures", "state", "all", "binCount", opt.binCountSeiz).avg;
        plotComparativeViolin({d1, d2}, xlabels);
        ylabel(sprintf('Avg seizures rate (Hz) - binsize : %d s', opt.binCountSeiz));
        title('Seizures'); grid on; set(gca, 'FontSize', 10);
        beautifyPlot();
        plotIdx = plotIdx + 1;

    else
        for j = 1:numel(state)
            subplot(nRows, nCols, plotIdx);
            ratsTM_ = ratTM{1}; ratsNPX_ = ratNPX{1};
            if state(j) ~= "all"
                ratsTM_ = ratTM{2};
                ratsNPX_ = ratNPX{2};
            end
            d1 = TM.getFrequency("rat", ratsTM_, "nature", "IS", "state", state(j), "binCount", opt.binCountIS).avg;
            d2 = NPX.getFrequency("rat", ratsNPX_, "nature", "IS", "state", state(j), "binCount", opt.binCountIS).avg;
            plotComparativeViolin({d1, d2}, xlabels);
            ylabel(sprintf('Avg IS rate (Hz) - binsize : %d s', opt.binCountIS));
            title(sprintf('IS (%s)', state(j))); grid on; set(gca, 'FontSize', 10);
            beautifyPlot();
            plotIdx = plotIdx + 1;
        end
    end
end

% Sleep vs wake comparison (pooled rats)
subplot(nRows, nCols, plotIdx)
dSleep = [TM.getFrequency("rat", ratTM{1}, "nature", "IS", "state", "sleep", "binCount", opt.binCountSeiz).avg; ...
          NPX.getFrequency("rat", ratNPX{1}, "nature", "IS", "state", "sleep", "binCount", opt.binCountSeiz).avg];
dWake  = [TM.getFrequency("rat", ratTM{1}, "nature", "IS", "state", "wake", "binCount", opt.binCountSeiz).avg; ...
          NPX.getFrequency("rat", ratNPX{1}, "nature", "IS", "state", "wake", "binCount", opt.binCountSeiz).avg];
plotComparativeViolin({dSleep, dWake}, {'Sleep','Wake'});
title('IS sleep vs wake - All rats pooled'); beautifyPlot();

% Global formatting
sgtitle("Distribution over our populations");

figName = "GlobalDistributionOverPopulations";
savePath = fullfile("/mnt/hubel-data-103/Corentin/Matlab/Figures_PSD", figName);
exportgraphics(fig, savePath + ".png", 'Resolution', 300);

end