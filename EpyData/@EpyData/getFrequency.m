function Count = getFrequency(this, opt)
% GETFREQUENCY - Compute time-binned frequency of events (IS or seizures) per rat 
%                depending on the state and compute the normalization time factors.
%
% Syntax:
%   Count = this.getFrequency(opt)
%
% Description:
%   Computes the frequency (Hz) of interictal spikes (IS) or seizures in fixed-size time bins.
%   Can restrict computation to a specific behavioral state ("sleep", "wake", or "all").
%   Frequencies are normalized based on time spent in the state per bin (only for IS).
%   If no events are found, an empty result is stored for the rat.
%
% Inputs:
%   this        - (1,1) EpyData instance
%   opt         - Options structure with fields:
%       .rat      : (1,:) double - Rat IDs to include (default: all)
%       .nature   : (1,1) string - "IS" or "seizures" (default: "IS")
%       .state    : (1,1) string - "all", "sleep", or "wake" (default: "all")
%       .binCount : (1,1) double - Time bin size in seconds (default: 3600)
%       .range    : (1,1) double - Max time (s) to consider (default: [] = full : go to the rangeCutIS)
%
% Output:
%   Count - Struct with fields:
%       .freq        : Frequency vectors (cell array, one per rat, NaNs if no data)
%       .avg         : Average frequency per rat (NaN-omitted)
%       .rat         : Rat IDs
%       .nature      : "IS" or "seizures"
%       .state       : "sleep", "wake", or "all"
%       .binSize     : Bin size in seconds
%       .normFactor  : Normalization factor applied (cell array per rat)
%
% Methods and functions used : 
%   - stateFractionPerBin
%   - IS

% Parse arguments
arguments
    this (1,1) EpyData
    opt.rat(:, 1) double = this.ratID
    opt.nature (1,1) string {mustBeMember(opt.nature, ["IS", "seizures"])} = "IS"
    opt.state (1,1) string {mustBeMember(opt.state, ["all", "sleep", "wake"])} = "all" % By default, all IS regardless of the state
    opt.binCount (1,1) double = 3600  % Bin size in seconds for the count
    opt.range (:, 1) double {mustBeNonnegative} = []
end

% Initialize output structure
Count = struct();
Count.nature = opt.nature;
Count.state = opt.state;
Count.binSize = opt.binCount;
Count.avg = 0;
idx = 0;

% Loop over rats to compute frequencies
for j = 1:numel(opt.rat)

    for i = 1:numel(this.ratID)
        if this.ratID(i) == opt.rat(j)
            idx = idx + 1;
            
            % Get event timestamps
            switch opt.nature
                case "seizures"
                    timestamps = this.seizuresTimestamps{i,1};

                    if isempty(timestamps)
                        fprintf("No seizures computed for the rat %d \n", opt.rat(j))
                        Count.freq{idx,1} = [];
                        Count.rat(idx,1) = opt.rat(j); % even if empty, we keep the mapping
                        break % Leave the for loop on i
                    end
                    timestamps = mean(timestamps, 2); % seizures are computed thanks to time intervals, to give a unique timestamp use mean of start/stop timestamps.

                case "IS"
                    timestamps = cell2mat(this.IS("rat", opt.rat(j), "state",opt.state));


            end

            if isempty(timestamps)
                fprintf("No sleep/wake computed for the rat %d \n", opt.rat(j))
                Count.freq{idx,1} = [];
                Count.rat(idx,1) = opt.rat(j); % even if empty, we keep the mapping
                break % Leave the for loop on i
            end

            % Define time range
            % Range of the count [0 range]. If no special windows, go until
            % the rangeCutIS.
            if  isempty(opt.range) || isempty(opt.range(j)) || opt.range(j) > max(timestamps)
                range = this.rangeCutIS(i);
                    % range = max(timestamps);
            else
                range = opt.range(j);
            end

            % Define bin edges and count events
            edges = 0:opt.binCount:floor(range)+opt.binCount;
            raw_counts = histcounts(timestamps', edges);

            % Compute the normalization factor regarding recorded time spent in state per bin
            normFactor = stateFractionPerBin(this,"rat", opt.rat(j), "nature", opt.nature, "state",opt.state,"range", range, "binCount", opt.binCount); % Special normalization factor depending the state

            % Adjust vector lengths if necessary
            if numel(normFactor) > 1 % We have a specific normalzation factor for every bin (sleep and wake)
                % Resize the data (at this stage, numel(normFactor) <= numel(raw_counts))
                raw_counts = raw_counts(:,1:numel(normFactor));
            end
            
            % Normalize and compute frequency
            freq = raw_counts' ./ normFactor; % Compute the count and convert it in frequency in Hz based on the time in the state per bin.
            % The vector freq have NaN element in case of no single seconds
            % in the state during a bin (can not have any event also in
            % that satate) so 0/0 even thought we have the knowledge that
            % there was indeed no timin that state
            % or also in case of no recording or no labeling of the state
            % (0/0)

            Count.freq{idx,1} = freq;
            Count.avg(idx,1) = mean(freq, 'omitnan'); % Average of frequency (normalized by binCount)
            Count.rat(idx,1) = opt.rat(j);
            Count.normFactor{idx,1} = normFactor; % Normalizations factor

            % Debug
            if any(freq == inf)
                warning("infinity -> need to be corrected")
                freq(freq == inf) = 0;
            end
            break
        end
    end
end

% Handle empty output
if idx == 0
    warning("No matching rat found in opt.rat. Returning empty output.");
    return;
end

end