function normFactor = stateFractionPerBin(this,opt) 
% Compute the time spent in a specific state per time bin for one rat.
%
% Syntax:
%   normFactor = stateFractionPerBin(this, opt)
%
% Description:
%   Calculates, for each time bin of specified duration, the total time (in seconds)
%   the rat spent in a specific state ("sleep" or "wake"). This is used to normalize
%   event frequencies in getFrequency.
%
% Inputs:
%   this        - (1,1) EpyData instance
%   opt         - Structure with options:
%     .rat      - (1,1) double : ID of the rat
%     .state    - (1,1) string : "all", "sleep", or "wake" (default: "all")
%     .binCount - (1,1) double : Bin size in seconds (default: 3600)
%     .range    - (1,1) double : Optional maximum time range (0 = use full duration) (default: 0)
%
% Output:
%   normFactor - Vector (nBins x 1) giving, for each bin, the amount of time
%                (in seconds) the rat spent in the specified state.
%
% Notes:
%   - If state = "all", returns a constant factor equal to bin size.
%   - Sleep/wake intervals are reshaped to avoid overlapping bins and ensure accurate bin assignment.
%   - A warning is raised if any bin has >100% occupancy (should not happen).


arguments
    this (1,1) EpyData
    opt.rat(1,1) double  % ID of the rat
    opt.nature (1,1) string {mustBeMember(opt.nature, ["IS", "seizures"])}
    opt.state (1,1) string {mustBeMember(opt.state, ["all", "sleep", "wake"])} = "all" % By default, all IS regardless of the state
    opt.binCount (1,1) double = 3600  % Bin size in seconds for the count
    opt.range (1,1) double {mustBeNonnegative} = 0
end

i = find(this.ratID == opt.rat);

% The time window of the entire recording
stateInt = [0, seconds(this.refDatesTimestamps(i,2) - this.refDatesTimestamps(i,1))];

switch opt.nature
    case "seizures"
        if size(this.unrecInt{i},1) > 0
            % construction of the recInt based on the unrecInt
            unrecInt = this.unrecInt{i};
            recInt = NaN(size(unrecInt,1)+1,2);
            window = [0, seconds(this.refDatesTimestamps(i,2) - this.refDatesTimestamps(i,1))];
            recInt(1,1) = window(1);
            recInt(2:end,1) = unrecInt(:,2);
            recInt(1:end-1,2) = unrecInt(:,1);
            recInt(end,2) = window(2);

            stateInt = recInt;
        end

    case "IS"
        switch opt.state
            case "all"
                if size(this.unrecInt{i},1) > 0
                    % construction of the recInt based on the unrecInt
                    unrecInt = this.unrecInt{i};
                    recInt = NaN(size(unrecInt,1)+1,2);
                    window = [0, seconds(this.refDatesTimestamps(i,2) - this.refDatesTimestamps(i,1))];
                    recInt(1,1) = window(1);
                    recInt(2:end,1) = unrecInt(:,2);
                    recInt(1:end-1,2) = unrecInt(:,1);
                    recInt(end,2) = window(2);

                    stateInt = recInt;
                end

            case {"wake","sleep"}
                if opt.state == "wake"
                    stateInt = this.wakeInt{i};
                else
                    stateInt = this.sleepInt{i};
                end
        end

end

%% Construction of the retaylored state intervals vector to intergrate bin size. Adding "artificial" edges (multiples of binCount) to cut overlapping intervals on two time bin.
binVector = repelem(opt.binCount*(1:ceil(stateInt(end,2)/opt.binCount)),2).'; % creation of a column vector containing all the multiple of Opt.binCount twice
flatVector = sortrows([stateInt(:); binVector]);
stateIntTayloredBinSize = [flatVector(1:2:end-1) , flatVector(2:2:end)];

%% Summing the duration of intervals of the state intervals vector per bin
bin_idx = floor(stateIntTayloredBinSize(:,1)./opt.binCount)+1; % Attributes an bin index for all the lines of the new state intervals vector (with no overlapping possible)
normFactor = accumarray(bin_idx, stateIntTayloredBinSize(:,2)-stateIntTayloredBinSize(:,1));

% Resize the normFactor vector based on the opt.range -> sometimes
% cut a lot (meaning the last timestamp of IS or seizure is really
% lower than the last timestamp of sleep). Sometimes too short
if ~opt.range == 0 && ceil(opt.range/opt.binCount) <= numel(normFactor)
    normFactor = normFactor(1:ceil(opt.range/opt.binCount),1);
end

%% Debugging Test
if any(normFactor>opt.binCount)
    error("IMPOSSIBLE NORM FACTOR, amount of state too important per bin")
end
end

        % % Indexes of state intervals overlapping two bins
        % startBin = floor(stateInt(:,1) / opt.binCount);
        % endBin = floor(stateInt(:,2) / opt.binCount);
        % idx = find(startBin ~= endBin);  % index where to split
        % 
        % % Creation of the small splitted intervals
        % cutTimes = opt.binCount * ceil(stateInt(idx,1) / opt.binCount);
        % firstPart = [stateInt(idx,1), cutTimes];
        % secondPart = [cutTimes, stateInt(idx,2)];
        % 
        % % Conservation of non overlapping intervals
        % mask = true(size(stateInt,1),1);
        % mask(idx) = false;
        % stateIntKept = stateInt(mask,:);
        % 
        % % Final construction gathering and sorting all the pieces intervalls
        % stateIntTayloredBinSize = sortrows([stateIntKept; firstPart; secondPart]);

        %% BACKUP

        % arguments
        % this (1,1) EpyData
        % opt.rat(1,1) double  % ID of the rat
        % opt.state (1,1) string {mustBeMember(opt.state, ["all", "sleep", "wake"])} = "all" % By default, all IS regardless of the state
        % opt.binCount (1,1) double = 3600  % Bin size in seconds for the count
        % opt.range (1,1) double {mustBeNonnegative} = 0
        % end
        % 
        % switch opt.state
        %     case "all"
        %         normFactor = opt.binCount; % No special normalization factor
        % 
        %     case {"wake","sleep"}
        %         i = find(this.ratID == opt.rat);
        % 
        %         if opt.state == "wake"
        %             stateInt = this.wakeInt{i};
        %         else
        %             stateInt = this.sleepInt{i};
        %         end
        % 
        %         %% Construction of the retaylored state intervals vector to intergrate bin size. Adding "artificial" edges (multiples of binCount) to cut overlapping intervals on two time bin.
        %         binVector = repelem(opt.binCount*(1:ceil(stateInt(end,2)/opt.binCount)),2).'; % creation of a column vector containing all the multiple of Opt.binCount twice
        %         flatVector = sortrows([stateInt(:); binVector]);
        %         stateIntTayloredBinSize = [flatVector(1:2:end-1) , flatVector(2:2:end)];
        % 
        %         %% Summing the duration of intervals of the state intervals vector per bin
        %         bin_idx = floor(stateIntTayloredBinSize(:,1)./opt.binCount)+1; % Attributes an bin index for all the lines of the new state intervals vector (with no overlapping possible)
        %         normFactor = accumarray(bin_idx, stateIntTayloredBinSize(:,2)-stateIntTayloredBinSize(:,1));
        % 
        %         % Resize the normFactor vector based on the opt.range -> sometimes
        %         % cut a lot (meaning the last timestamp of IS or seizure is really
        %         % ower than the last timestamp of sleep). Sometimes too short
        %         if ~opt.range == 0 && ceil(opt.range/opt.binCount) <= numel(normFactor)
        %             normFactor = normFactor(1:ceil(opt.range/opt.binCount),1);
        %         end
        % 
        %         %% Debugging Test
        %         if any(normFactor>opt.binCount)
        %             warning("IMPOSSIBLE NORM FACTOR, amount of state too important per bin")
        %         end
        % end
        % 
        % end

