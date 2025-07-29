function plotStatPowerSpectralDensity(this, opt)
% plotStatPowerSpectralDensity - Plot and test the statistical power spectrum of IS occurrences.
%
% This method compares the power spectral density (PSD) of interictal spike (IS) occurrences
% with surrogate data to identify frequency bands (e.g., multidian) where the real data shows
% significantly higher power. The surrogate data preserve the circadian structure via block shuffling.
%
% INPUTS:
%   opt.rat        - (1,:) double. Rat ID(s) to include. Default = all rats in object.
%   opt.state      - (1,1) string. Behavioral state ("all", "sleep", "wake"). Default = "all".
%   opt.binCount   - (1,1) double. Time bin size (s) used for IS histogram. Default = 3600.
%   opt.range      - (1,1) double. Time range for the analysis (days). Default = 0 (use full session).
%   opt.nSurro     - (1,1) double. Number of surrogates. Default = 100.
%   opt.binShuffle - (1,1) double. Shuffling block size (s). Default = 86400 (1 day).
%
% This function also updates the property `fMultidianCycl` with the dominant frequency
% (in Hz) if existing of significant multidian peaks for each rat and state.


arguments
    this (1,1) EpyData
    opt.rat(1,:) double = this.ratID  % optionnal, by défault all the rats
    opt.state (1,1) string {mustBeMember(opt.state, ["all", "sleep", "wake"])} = "all" % By default, all IS regardless of the state
    opt.binCount (1,1) double = 3600  % Bin size in seconds for the count of the data of interests
    opt.range (1,1) double {mustBeNonnegative} = 0
    opt.nSurro (1,1) double {mustBeNonnegative, mustBeInteger} = 100 % Number of surrogates expected by rats If 0, no surrogate, work with the real data.
    opt.binShuffle (1,1) double {mustBePositive} = 3600*24 % Bin size in seconds used to shuffle blocks of data count to create the surrogates. Default 1 day to keep circadian cyclicity.
end


%% Data
PSDReal = this.computePSD("rat", opt.rat, "state", opt.state, "range", opt.range, "binCount",opt.binCount, "nSurro",0);
PSDSurro = this.computePSD("rat", opt.rat, "state", opt.state, "range", opt.range ,"binCount", opt.binCount, "nSurro",opt.nSurro, "binShuffle", opt.binShuffle);

sigInt = cell(numel(PSDReal.pxx),1);

% Boundaries for statistics 
freqBounds = [1/(24*12*3600) 1/(3*24*3600)]; % Statistics beetween 3 and 12 days

for i=1:numel(PSDReal.pxx)
 
    idxReal = PSDReal.f{i} >= freqBounds(1) & PSDReal.f{i} <= freqBounds(2);
    idxSurro = PSDSurro.f{i} >= freqBounds(1) & PSDSurro.f{i} <= freqBounds(2);

    fSize = min(numel(PSDReal.f{i}(idxReal)),numel(PSDSurro.f{i}(idxSurro)));
    numBinStat = floor(fSize/3); % 3 points of PSD per bin to do the stat
    binEdges = linspace(freqBounds(1), freqBounds(2), numBinStat + 1);
    binIndicesReal = arrayfun(@(lo, hi) find(PSDReal.f{i} >= lo & PSDReal.f{i} < hi), binEdges(1:end-1), binEdges(2:end), 'UniformOutput', false); % Get frequency indices for each bin as a cell array 
    binIndicesSurro = arrayfun(@(lo, hi) find(PSDSurro.f{i} >= lo & PSDSurro.f{i} < hi), binEdges(1:end-1), binEdges(2:end), 'UniformOutput', false); % Get frequency indices for each bin as a cell array 

    position_real = zeros(numBinStat, 1); % -1 (real significatively below surrogate), 0 (not significative), 1 (real significatively above surrogate)

    for j = 1:numBinStat
        iReal = binIndicesReal{j}; % Extract indices of frequencies belonging to bin i
        iSurro = binIndicesSurro{j};

        sorted_surro = sort(mean(PSDSurro.pxx{i}(iSurro, :), 1)); % Mean vector in croissant order of surrogate PSD in bin. Each column is the mean of 1 surrogate.

        real_mean_bin = mean(PSDReal.pxx{i}(iReal)); % Mean value of real data in bin

        % Get empirical rank of the real data if it would have been put inside
        % the surrogate mean vector
        rank_real = find(sorted_surro >= real_mean_bin, 1);
        if isempty(rank_real)
            rank_real = length(sorted_surro) + 1;
        end

        p_val = rank_real / (length(sorted_surro) + 1); % Probability of the position. If 0, real data below all the surrogate mean vector values, if 1 above all of them.

        % Double-tailed test at alpha = 0.05
        if p_val <= 0.025 % The first quartile
            position_real(j) = -1;
        elseif p_val >= 0.975 % The last quartile
            position_real(j) = 1;
        end

    end
    sigInt{i} = position_real;


%% Plot

figure(Position=get(0,"ScreenSize"));
f_max = 1/opt.binCount;
interestfReal = PSDReal.f{i} <= f_max;
plot(log(1./(PSDReal.f{i}(interestfReal)*3600*24)), (PSDReal.pxx{i}(interestfReal))) % fréquence en cycles/jour (ou autre unité)
hold on
interestfSurro = PSDSurro.f{i} <= f_max;
semplot(log(1./(PSDSurro.f{i}(interestfSurro,:)*3600*24)), (PSDSurro.pxx{i}(interestfSurro,:)))
xlim([-1, log(15)]);
xticks([log(1/2), log(1), log(3), log(7), log(10), log(15)]);
xticklabels({'12h', '24h', '3d', '7d', '10d', '15d'});
title(sprintf("Real vs Surrogate PSD - %s IS – Rat %d –", opt.state, opt.rat(i)));
xlabel('Time');
ylabel('Power/Frequency (AU/days⁻¹)');
grid on;
hold on;
% Plot significance bars for bins where real PSD > surrogate
y_max = max(PSDReal.pxx{i}(interestfReal)) * 1.1;

for j = 1:numBinStat
    if sigInt{i}(j) == 1
        % Convert bin edge frequencies to log time units
        f_high = binEdges(j + 1);
        f_low = binEdges(j);
        x_high = log(1 / (24 * 3600 * f_low));
        x_low = log(1 / (24 * 3600 * f_high));

        plot([x_low, x_high], [y_max, y_max], '-k', 'LineWidth', 3);
        hold on;
    end
end

%% Seak of the multidian pic for all IS 
% Find index of this rat in the class property
idx = find(this.ratID == opt.rat(i));
if isempty(idx)
    warning("Rat %d not found in the object's ratID list.", opt.rat(i));
    continue;
end

% Initialize with NaN in case no significant bins are found
peakFreq = NaN;

% Identify bins where real PSD is significantly higher than surrogate
valid_bins = find(sigInt{i} == 1);

if ~isempty(valid_bins)
    freq_idxs = [];
    for j = valid_bins'
        % Collect frequency indices from all significant bins
        freq_idxs = [freq_idxs, binIndicesReal{j}];
    end
    freq_idxs = unique(freq_idxs);  % Remove duplicates if any

    % Find the frequency index with the maximum PSD value in those bins
    [~, idx_max] = max(PSDReal.pxx{i}(freq_idxs));
    peakFreq = PSDReal.f{i}(freq_idxs(idx_max));  % Peak frequency
end

% Store the result in the class property for this rat
switch opt.state
    case "all"
        this.fMultidianCycl(idx, 1) = peakFreq;
    case "sleep"
        this.fMultidianCycl(idx, 2) = peakFreq;
    case "wake"
        this.fMultidianCycl(idx, 3) = peakFreq;
end

end
end

