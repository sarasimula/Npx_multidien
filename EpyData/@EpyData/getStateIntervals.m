function stateInt = getStateIntervals(this, opt)
% GETSTATEINTERVALS - Returns time intervals corresponding to a given brain state for selected rats.
%
% Syntax:
%   stateInt = getStateIntervals(this)
%   stateInt = getStateIntervals(this, opt)
%
% Description:
%   Retrieves the intervals corresponding to the specified brain state ("sleep", "wake", "unrec", or "unknown")
%   for each specified rat. The "unknown" state is defined as the complement of all known intervals within the
%   recording window.
%
% Inputs:
%   this        - (EpyData) Object containing state interval data and reference timestamps.
%   opt         - (struct, optional) Name-value options:
%       • rat     : [1×N double] Rat IDs to include (default = all rats in object)
%       • state   : ["sleep"|"wake"|"unrec"|"unknown"] Brain state of interest (default = "sleep")
%
% Outputs:
%   stateInt    - Cell array [nRats×1], each cell contains a [M×2] double matrix of intervals [start, end] in seconds
%                 since session start, corresponding to the selected state for each rat.
%
% Notes:
%   - For "unknown", the function returns the intervals not covered by any of the three known states
%     (sleep, wake, unrec) within the total recording duration.
%   - If a specified rat is not found, a warning is issued and the corresponding cell is left empty.
%
% See also: refDatesTimestamps, sleepInt, wakeInt, unrecInt

arguments
    this (1,1) EpyData
    opt.rat(1,:) double = this.ratID  % optionnal, by défault all the rats
    opt.state (1,1) string {mustBeMember(opt.state, ["sleep", "wake", "unknown", "unrec"])} = "sleep" 
end

stateInt = cell(numel(opt.rat),1);   % Initialize output cell array
idx = 0;                            % Index to track output cell position

% Loop over each requested rat ID
for j = 1:numel(opt.rat)
    previousIdx = idx;  % Store previous index to detect if rat was found
    
    % Loop over all rats in the object to find matching rat ID
    for i = 1:numel(this.ratID)
        if this.ratID(i) == opt.rat(j)
            idx = idx + 1;  % Increment output index
            
            % Select intervals depending on the requested brain state
            switch opt.state
                case "sleep"
                    stateInt{idx} = this.sleepInt{i};  % Sleep intervals for this rat
                case "wake"
                    stateInt{idx} = this.wakeInt{i};   % Wake intervals for this rat
                case "unrec"
                    stateInt{idx} = this.unrecInt{i};  % Unrecognized intervals for this rat
                case "unknown"
                    % Concatenate all known intervals and sort by start time
                    allInt = sortrows([this.sleepInt{i}; this.wakeInt{i}; this.unrecInt{i}]);
                    
                    % Define the full recording time window for this rat
                    windowTime = [0, seconds(this.refDatesTimestamps(i,2) - this.refDatesTimestamps(i,1))];
                    
                    % Initialize matrix for unknown intervals (one more row than known intervals)
                    unknownInt = zeros(size(allInt,1)+1, 2);
                    unknownInt(1,1) = windowTime(1);           % Start at beginning of session
                    unknownInt(end,2) = windowTime(2);         % End at end of session
                    
                    % Construct unknown intervals as gaps between known intervals
                    unknownInt(2:end,1) = allInt(:,2);         % Start times are ends of known intervals
                    unknownInt(1:end-1,2) = allInt(:,1);       % End times are starts of known intervals
                    
                    % Filter intervals with duration > 1 nanosecond (avoid zero-length)
                    validIdx = find(unknownInt(:,2) - unknownInt(:,1) > 1e-9);
                    stateInt{idx} = unknownInt(validIdx, :);
            end
        end
    end
    
    % If no rat found, issue a warning and leave output cell empty
    if previousIdx == idx
        warning("No matching rat found for rat ID %d. Returning empty output.", opt.rat(j));
    end
end
end