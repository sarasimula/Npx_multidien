function this = analyzePSD(this, opt)
% ANALYZEPSD - Analyze Power Spectral Density and extract multidien features
%
% This method estimates the Power Spectral Density (PSD) of interictal spike (IS) rates
% across animals and vigilance states, compares it to surrogate data, identifies
% frequency bands in the multidien range (3–15 days) with significant power increases,
% and extracts peak frequencies and associated metrics.
%
% Updated object properties:
%   - this.PSDsurro : Cell array storing surrogate PSDs for each rat and vigilance state
%   - this.multidienCycl{:,1}.f0 : (1×n double) Significant peak frequencies [Hz] in the multidian range
%                                  (3 to 15 days, or up to the duration of the session) 
%                                  where real PSD is significantly higher than surrogates 
%                                  (computed using the 'all' state).
%                                  Ordered in ascendant order of importance of the peak.
%   - this.multidienCycl{idx,1}.importanceFact_f : (Vector 1xn) Relative importance of the peak associated to the frequency of f0
%                                   (the more it's close to 1 the more it's important).
%                                   Ordered in the same order as f0.
%   - this.multidienCycl{idx,1}.fWidth : (Vector nx2)  Boundaries frequency [fStart, fEnd] describing each peak at f0. Calculated finding the f associated to power_f0/2.
%   - this.multidienCycl{idx,1}.sigBandEdges : (Vector mx2)  Frequency bins [fStart, fEnd] significantly above surrogates.
%   - this.multidienCycl{idx,1}.freqISFiltered : structure with field of all state (.all, .sleep, .wake) giving the multidien temporal serie rhythm (after band passing at each f0).
%   - this.multidienCycl{:,1}.incrFromShufSleep : (Vector 1xn) Each increase associated to a f0 peak, between normalized real PSD and 
%                                  surrogate mean in the range of the peak width, in the 'sleep' state
%   - this.multidienCycl{:,1}.incrFromShufWake : Same as above, for the 'wake' state
%   - this.multidienCycl{:,1}.areaFromShufSleep : (Vector 1xn) Each area associated to a f0 peak, between normalized real PSD and 
%                                  surrogate mean in the range of the peak width, in the 'sleep' state
%   - this.multidienCycl{:,1}.areaFromShufWake : Same as above, for the 'wake' state
%
% INPUTS (via `opt`):
%   - rat        : (1,:) double - IDs of rats to analyze (default = all rats)
%   - state      : string - 'all' | 'sleep' | 'wake' (default = 'all')
%   - binCount   : double - Bin size (s) for IS histogram used to compute PSD (default = 3600)
%   - range      : double - Analysis duration (days), 0 = full session (default = 0)
%   - nSurro     : integer - Number of surrogate datasets to generate (default = 100)
%   - binShuffle : double - Bin size (s) used to shuffle surrogate data (default = 86400)
%   - force      : logical - Force recomputation of PSD even if it exists (default = false)
%
%
% OUTPUTS:
%   - this         : Updated object with `fMultidianCycl` set for each rat and state
%
% Methods and functions used : 
%   - computePSD
%   - findpeaks_Matlab, the built-in function findpeaks of matlab and not
%   the one from the FMAT toolbox (conflict)

% Parse inputs
arguments
this (1,1) EpyData
opt.rat (1,:) double = this.ratID
opt.state (1,1) string {mustBeMember(opt.state, ["all", "sleep", "wake"])} = "all"
opt.binCount (1,1) double = 3600
end

% Initialization output

% Frequency range to do the stats corresponding to 3–15 day cycles
% (multidian). The range is actualized after to take in account time of rec.
switch this.samplingFreq
    case 512 % TM
        freqBounds = [1/(21*24*3600), 1/(3*24*3600)]; % [f21 days, f3 days]                
    case 30000 % NPX
        freqBounds = [1/(15*24*3600), 1/(3*24*3600)]; % [f15 days, f3 days]
end

% Loop over rats
for i = 1:numel(opt.rat)
    ratID = opt.rat(i);
    idx = find(this.ratID == ratID);
    if isempty(idx)
        warning("Rat %d not found in this.ratID", ratID);
        continue;
    end

    % Access real and surrogate PSDs
    pxxReal = this.PSDreal{idx}.(opt.state).pxx;
    fReal   = this.PSDreal{idx}.(opt.state).f;
    pxxSurro = this.PSDsurro{idx}.(opt.state).pxx;
    fSurro   = this.PSDsurro{idx}.(opt.state).f;

    %% Sanity check
    if fReal ~= fSurro
        warning("Be careful with the two frequencies grid")
    end
    minfOfTheSignal = 1/seconds(this.refDatesTimestamps(idx,2)- this.refDatesTimestamps(idx,1)); % minimum frequency can be given by the time range of the recording if very short
    freqBounds(1) = max(freqBounds(1),minfOfTheSignal);               

    % Normally well done upper in the pipeline : but select acceptable frequencies in the multidian range
    idxReal = fReal >= freqBounds(1) & fReal <= freqBounds(2);
    idxSurro = fSurro >= freqBounds(1) & fSurro <= freqBounds(2);

    %%
    % Split into frequency bins for statistical testing
    nFreqs = min(sum(idxReal), sum(idxSurro));
    nBins = max(floor(nFreqs / 3)); % 3 points of PSD per bin to do the stat
    binEdges = linspace(freqBounds(1), freqBounds(2), nBins + 1);

    % Initialize significance vector
    sigVec = zeros(nBins,1); % -1 (real significatively below surrogate), 0 (not significative), 1 (real significatively above surrogate)
    binReal = cell(nBins,1);
    binSurro = cell(nBins,1);

    for b = 1:nBins
        % Extract indices of frequencies belonging to bin i
        binReal{b} = find(fReal >= binEdges(b) & fReal < binEdges(b+1));
        binSurro{b} = find(fSurro >= binEdges(b) & fSurro < binEdges(b+1));

        if isempty(binReal{b}) || isempty(binSurro{b})
            continue
        end

        % Compute mean power in bin
        meanReal = mean(pxxReal(binReal{b})); % Mean value of real data in bin
        meanSurro = sort(cellfun(@(p) mean(p(binSurro{b})), pxxSurro)); % Mean vector in croissant order of surrogate PSD in bin. Each column is the mean of 1 surrogate.

        % Empirical p-value
        if sum(meanReal >= meanSurro) > 0
            rank = sum(meanReal >= meanSurro) + 1;
        else
            rank = 0;
        end
        
        p = rank / (numel(meanSurro) + 1); % Probability of the position. If 0, real data below all the surrogate mean vector values, if 1 above all of them.

        % Double-tailed test at alpha = 0.05
        if p >= 0.975
            sigVec(b) = 1; % The first quantile % Real PSD significantly above surrogates
        elseif p <= 0.025
            sigVec(b) = -1; % The last quantile % Real PSD significantly below surrogates
        end
    end

    % Identify peaks in significant positive bins
    [peakFreq,peakValsSorted] = deal([]);
    sigBins = find(sigVec == 1);


        %%              %         Compute and update the fMultidianCycl property of EpyData class based only on state = "all"        %
      
  if opt.state == "all"
        % Reset the property fMultidianCycl
        this.multidienCycl{idx,1}.f0 = [];
        this.multidienCycl{idx,1}.importanceFact_f = [];
        this.multidienCycl{idx,1}.fWidth = [];
    
        % Compute the frequencies of multidien rythm
        if ~isempty(sigBins)

            % Output of the function 
            sigStart = binEdges(sigBins);
            sigEnd = binEdges(sigBins+1);
            this.multidienCycl{idx,1}.sigBandEdges = [sigStart(:), sigEnd(:)];

            % Locate peaks within significant bins using findpeaks
            idxSig = unique(cell2mat(binReal(sigBins)));
            fSig = fReal(idxSig);
            pxxSig = pxxReal(idxSig);

            [peakVals, locs] = findpeaks_Matlab(pxxSig, fSig);

            % Order the pics in descendant values
            [peakValsSorted, sortIdx] = sort(peakVals, 'descend');
            peakFreq = locs(sortIdx);                          % locs ordered based on peakVals
            peakValsSorted = peakValsSorted / max(pxxReal);    % Normalization

                    % % Former method when only one peak :  Compute the new max peak
                    % if ~isempty(sigBins)
                    %     % Record frequency edges of significant bins
                    %     sigStart = binEdges(sigBins);
                    %     sigEnd = binEdges(sigBins+1);
                    %     sigBandEdges{idx,1} = [sigStart(:), sigEnd(:)];
                    % 
                    %     % Locate max power frequency within significant bins
                    %     idxMax = unique(cell2mat(binReal(sigBins)));
                    %     [~, iMax] = max(pxxReal(idxMax));
                    %     peakFreq = fReal(idxMax(iMax));
                    % end
        end

        % Store in fMultidianCycl
        this.multidienCycl{idx,1}.f0 = peakFreq;
        this.multidienCycl{idx,1}.importanceFact_f = peakValsSorted;

        % Compute the width of peaks for band pass of filter
        if ~isempty(peakFreq)
             for F = 1:numel(peakFreq)
                 f0 = peakFreq(F);
                 [~, f0Idx] = min(abs(fReal - f0)); % Index of the pic
                 powerHalf = pxxReal(f0Idx) / 2; % Half of the peak amplitude

                 % Search the left cut frequency
                 i1 = f0Idx;
                 while i1 > 1 && pxxReal(i1) > powerHalf
                     i1 = i1 - 1;
                 end

                 % Search the right cut frequency
                 i2 = f0Idx;
                 while i2 < numel(pxxReal) && pxxReal(i2) > powerHalf
                     i2 = i2 + 1;
                 end

                 this.multidienCycl{idx,1}.fWidth(F,:) = [fReal(i1), fReal(i2)];
             end  
        end
  end

    %%              %         Compute the best slow-oscillation signal based on the result of the PSD         %
    fMultidien = this.multidienCycl{idx,1}.f0;
    if ~isempty(fMultidien)

        freq = this.getFrequency("rat", ratID, "nature", "IS", ...
            "state", opt.state, "binCount", opt.binCount).freq{1}';

        freqProcessed = signalProcessing(freq, nHarmonics=20); % NaN interpolation necessary to apply band pass filter 
            % (signalProcessing applies also smooth, eventhough not necessary, but without effect) 

        filters = cell(numel(fMultidien),1);
        fs = 1/ opt.binCount;
        for F = 1:numel(fMultidien)
            fWidth = this.multidienCycl{idx,1}.fWidth(F,:);
            filters{F} = designfilt('bandpassiir', 'FilterOrder', 4, ...
                'HalfPowerFrequency1', fWidth(1), ...
                'HalfPowerFrequency2', fWidth(2), ...
                'SampleRate', fs);
        end

        % Applicate multi-band filtering
        freqFiltered = zeros(size(freqProcessed));
        for k = 1:numel(filters) 
            freqFiltered = freqFiltered + filtfilt(filters{k}, freqProcessed);
        end

        % Normalization
        scaleFactor = std(freqProcessed) / std(freqFiltered);
        if isfinite(scaleFactor) && scaleFactor > 0
            freqFiltered = freqFiltered * scaleFactor + mean(freqProcessed)- mean(freqFiltered);
        end

        this.multidienCycl{idx,1}.freqISFiltered.(opt.state) = freqFiltered;

    else
        this.multidienCycl{idx,1}.freqISFiltered = [];

    end

                                % % Former method when only one peak :  Compute the new max peak
                                % if opt.state == "all"
                                %     freq = this.getFrequency("rat", ratID, "nature", "IS", ...
                                %         "state", opt.state, "binCount", opt.binCount).freq{1}'; % * opt.binCount;
                                %     t = (0:numel(freq)-1) / 24; % time in days
                                %     freqFiltered = signalProcessing(freq, nHarmonics=20);
                                %     f0 = this.multidienCycl{idx,1}.f0;
                                %
                                %     if ~isnan(f0)
                                %         omega = 2 * pi * f0 * (24*3600); % in d⁻¹
                                %
                                %         % Model : A*cos(ωt + φ) = a*cos(ωt) + b*sin(ωt)
                                %         X = [cos(omega*t(:)) sin(omega*t(:))];
                                %         params = X \ freqFiltered(:);
                                %
                                %         % Store results
                                %         this.multidienCycl{idx,1}.A = hypot(params(1), params(2));        % amplitude
                                %         this.multidienCycl{idx,1}.phi_t0 = atan2(-params(2), params(1));     % phase
                                %     else
                                %         this.multidienCycl{idx,1}.A = NaN;
                                %         this.multidienCycl{idx,1}.phi_t0 = NaN;
                                %     end
                                % end

    %%              %          Give the increase from shuffle of the sleep and wake data          %                   %
    increaseFromShuffle = [];
    areaFromShuffle = [];
    increaseFromShuffleString = "incrFromShuf" + opt.state.extractBetween(1,1).upper + extractAfter(opt.state,1);
    areaFromShuffleString = "areaFromShuf" + opt.state.extractBetween(1,1).upper + extractAfter(opt.state,1);

    % New method
    if opt.state ~= "all" && ~isempty(fMultidien)  
        for F = 1:numel(fMultidien)
            % Défine the frequence band based on the pic width
            fStart = this.multidienCycl{idx,1}.fWidth(F, 1);
            fEnd = this.multidienCycl{idx,1}.fWidth(F, 2);

            % Index of peak frequency and frequencies in the band of the peak
            idxPeakFreq = find(fReal == fMultidien(F));
            idxBandFreq = find(fReal >= fStart & fReal <= fEnd);

            % Extract surrogate pxx and real pxx in the peak and compute increaseFromShuffle at the peak
            meanPeakFreqNormSurro = mean(cellfun(@(p) p(idxPeakFreq)./max(p),pxxSurro)); % mean of the normalized power of surrogates at the peak frequency
            increaseFromShufflePeak = pxxReal(idxPeakFreq)/max(pxxReal) - meanPeakFreqNormSurro; % compute the diff beetween reel noramlized value and normalized mean surrogates.

            % Extract surrogate pxx and real pxx in the band and compute cumulative increaseFromShuffle in the band to compute area  
            increaseFromShuffleBand = 0;
            for idxF = 1:numel(idxBandFreq)
                meanAreaNormSurro = mean(cell2mat(cellfun(@(p) p(idxBandFreq(idxF))./max(p),pxxSurro, 'UniformOutput', false))); % mean of the normalized power area of surrogates around peak frequency
                increaseFromShuffleBand = increaseFromShuffleBand + pxxReal(idxBandFreq(idxF))/max(pxxReal) - meanAreaNormSurro; % Area = sum
            end
            
            increaseFromShuffle = [increaseFromShuffle, increaseFromShufflePeak];
            areaFromShuffle = [areaFromShuffle, increaseFromShuffleBand];
        end

    this.multidienCycl{idx,1}.(areaFromShuffleString) = areaFromShuffle;
    this.multidienCycl{idx,1}.(increaseFromShuffleString) = increaseFromShuffle;
    
    elseif opt.state ~= "all" && isempty(fMultidien)
        this.multidienCycl{idx,1}.(areaFromShuffleString) = areaFromShuffle;
        this.multidienCycl{idx,1}.(increaseFromShuffleString) = increaseFromShuffle;
    end
end


                                % % Former method when evaluation on an unique frequency
                                % if ~isnan(this.multidienCycl{idx,1}.f0)
                                %     peakFreq = this.multidienCycl{idx,1}.f0;
                                %     idxPeakFreq = find(fReal == peakFreq);
                                %     meanPeakFreqNormSurro = mean(cellfun(@(p) p(idxPeakFreq)./max(p),pxxSurro)); % mean of the normalized power of surrogates at the peak frequency
                                %     increaseFromShuffle = pxxReal(idxPeakFreq)/max(pxxReal) - meanPeakFreqNormSurro; % compute the diff beetween reel noramlized value and normalized mean surrogates.
                                % else, increaseFromShuffle = NaN;

