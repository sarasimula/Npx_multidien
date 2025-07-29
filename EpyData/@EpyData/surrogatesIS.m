function surro = surrogatesIS(this, opt)
%SURROGATESIS Generate surrogate interictal spike (IS) frequency time series.
%
%   surro = SURROGATESIS(obj) returns the true IS frequency time series for
%   all rats and all states, binned in 1-hour intervals by default.
%
%   surro = SURROGATESIS(obj, opt) allows specification of several options
%   for rat selection, brain state, bin size, and surrogate generation.
%   If surrogate generation is requested, the function returns a shuffled
%   version of the IS time series by permuting large time blocks, preserving
%   long-term structure such as circadian rhythms.
%
%   This method is particularly useful to assess the statistical relevance
%   of observed IS fluctuations by comparison with null distributions that
%   preserve coarse temporal dynamics but destroy fine-grained structure.
%
%   INPUT:
%     obj         : EpyData object containing telemetry or Neuropixels data.
%     opt         : (optional) struct with the following fields:
%       - rat         : [N x 1] double array of rat IDs. Default: all available.
%       - state       : string, "all" | "sleep" | "wake". Default: "all".
%       - binCount    : scalar, bin size in seconds to count the IS frequency. Default: 3600 (1 hour).
%       - range       : scalar, total duration in seconds to extract ([] = full : go to rangeCutIS). Default: [].
%       - nSurro      : integer ≥ 0, number of surrogates to generate per rat.
%                       If 0, return real IS frequency only. Default: 100.
%       - binShuffle  : scalar > 0, size in seconds of blocks to permute when
%                       generating surrogates. Default: 86400 (1 day).
%
%   OUTPUT:
%     surro      : Cell array of size [N x 1], one cell per rat.
%                  If nSurro > 0: each cell contains a [T x nSurro] matrix,
%                  where T is the number of bins.
%                  If nSurro == 0: each cell contains a [T x 1] vector of the real IS rate.
%
%   Dependencies : getFrequency method

arguments
    this (1,1) EpyData
    opt.rat(1,:) double = this.ratID  % By default all the rats
    opt.state (1,1) string {mustBeMember(opt.state, ["all", "sleep", "wake"])} = "all" % By default, all IS regardless of the state
    opt.binCount (1,1) double = 3600  % Bin size in seconds for the count of the data of interests
    opt.range (:, 1) double {mustBeNonnegative} = [] % Go to rangeCutIS
    opt.nSurro (1,1) double {mustBeNonnegative, mustBeInteger} = 100 % Number of surrogates expected by rats. If 0, no surrogate, just return the real ISfrequency.
    opt.binShuffle (1,1) double {mustBePositive} = 3600*24 % Bin size in seconds used to shuffle blocks of data count to create the surrogates. Default 1 day to keep circadian cyclicity.
end

% Preallocate output cell array (one entry per rat)
surro = cell(numel(opt.rat),1);

% Compute real IS frequency time series (already normalized by normFactor), binned over time
ISfrequency = this.getFrequency("rat",opt.rat,"nature","IS","state", opt.state, "binCount",opt.binCount, "range",opt.range).freq;

% Possibility to process (interpolation + smoothing) or no the signal before comparing 
% it to the surrogates. Methodology kept : we keep the real count without processing
ISprocessed = ISfrequency;
    % ISprocessed = cellfun(@(x) signalProcessing(x),ISfrequency, 'UniformOutput',false);

if opt.nSurro > 0
    % Loop over each rat to generate surrogates
    for i = 1:numel(ISprocessed)
        vecToShuffle = ISprocessed{i};
        n = numel(vecToShuffle);
        binPerShuffle = round(opt.binShuffle / opt.binCount);

        % Compute how many full blocks can be shuffled
        nBins = floor(n / binPerShuffle);
        if nBins == 0
            warning("Not enough data to shuffle for rat %d", opt.rat(i));
            continue
        end

        % Reshape the time series into columns of contiguous time blocks
        reshaped = reshape(vecToShuffle(1:nBins*binPerShuffle), binPerShuffle, nBins);
        nLeftover = n - nBins * binPerShuffle;

        % Generate surrogate permutations (random order of blocks)
        perms = zeros(nBins, opt.nSurro);
        for s = 1:opt.nSurro
            perms(:,s) = randperm(nBins)';
        end

        % Reconstruct surrogate time series based on permuted blocks
        surrogates = zeros(opt.nSurro, n);
        for s = 1:opt.nSurro
            shuffledBins = reshaped(:, perms(:,s));
            reshapedVec = shuffledBins(:);

            % Append leftover samples to a random position
            if nLeftover > 0
                r = randi(nBins);
                leftover = reshaped(1:nLeftover, r);
                reshapedVec = [reshapedVec; leftover];
            end

            surrogates(s,:) = reshapedVec;
        end

        % Store result: one matrix [T x nSurro] per rat
        surro{i,:} = surrogates';
    end

else
    % If no surrogate requested (opt.nSurro = 0), just return the processed real IS frequency 
    surro = ISprocessed;
end

end



