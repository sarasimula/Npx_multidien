function plotMultidienDistr(this, opt)
% PLOTMULTIDIENDISTR — Plot the distribution of multidien cycle periods across rats
%
% This method visualizes the distribution of the detected multidien cycles (in days)
% across a subset of rats, for a given vigilance state (all, sleep, or wake).
%
% INPUTS:
%   this        — an instance of the EpyData class
%   opt.rat     — vector of rat IDs to include in the plot (default: all in this.ratID)
%
% The function will ignore values equal to 0 (no computation) or NaN (no significant cycle),
% and will issue warnings accordingly. It automatically infers the data modality (TM or NPX)
% based on the sampling frequency.


arguments
    this (1,1) EpyData
    opt.rat(1,:) double = this.ratID
end

% Determine which type of instance
switch this.samplingFreq
    case 512 % TM
        modality = "TM";
    case 30000 % NPX
        modality = "NPX";
end

% project only on the ratID of interest
idx = [];
for i = 1:numel(opt.rat) % Select indices corresponding to the specified rats
    ratID = opt.rat(i);
    idx = [idx, find(this.ratID == ratID)];
    if isempty(idx)
        warning("Rat %d not found in this.ratID", ratID);
        continue;
    end
end
multidienCyclOfInterest = [this.multidienCycl(idx)]; 
noEmptyMultidienCycl = multidienCyclOfInterest(cellfun( @(x) ~isempty(x.f0), multidienCyclOfInterest));
% Extract the main frequency of the multidien periodicity in the state "all"
fMultidien = cellfun (@(x) x.f0(1), noEmptyMultidienCycl);

% Convert frequencies to periods in days and define histogram bins
minDay = floor(min(1./fMultidien/(24*3600)));
maxDay = ceil(max(1./fMultidien/(24*3600)));
edges = minDay : maxDay;

% Plot settings
figure('Units', 'pixels', 'Position', [0, 0, 2400, 1600], ...
    'Visible','on');
histogram(1./fMultidien/(24*3600), BinEdges=edges);
title(sprintf('%s - Distribution of the period associated to the main multidien peak', modality));
xlabel('Time period (in days)');
ylabel('Number of rats');
ylims = ylim;
yticks(floor(ylims(1)) : ceil(ylims(2)));
xticks(minDay : maxDay);


% --- Save figure ---
% figName = "TM_histogramMultidienPeak";
% savePath = fullfile("/mnt/hubel-data-103/Corentin/Matlab/Figures_PSD", figName);
% exportgraphics(gcf, savePath + ".png", 'Resolution', 300);

end