function this = getCellAssemblies(this, opt)
%GETCELLASSEMBLIES Detects neuronal assemblies using ISAC method
%
%   This method identifies cell assemblies from spike data using the
%   ISAC approach. The detection is done separately for each rat and
%   session over a specified and checked global time window. The resulting
%   assemblies are stored in the `compositionISAC` field under the `assemblies`
%   property of the EpyData object.
%
%   Syntax:
%     this = getCellAssemblies(this, opt)
%
%   Inputs:
%     this                           - (1x1 EpyData) Instance of EpyData class, must contain
%                                      a non-empty `unitsAll` property.
%
%     opt (Name-Value Arguments):
%       .rat                         - (1xN double) IDs of rats to process (default: all rats in `this.ratID`)
%       .timeWindow2DetectAssemb     - (1x1 double) Size of the sliding time window [in seconds] for co-activation 
%                                      detection using ISAC (default: 0.03 s = 30 ms)
%
%   Changed Properties:
%     this.assemblies{r}.compositionISAC{s} 
%       - Cell array of length equal to the number of sessions. Each cell contains a 
%         1xA cell array (A = number of assemblies), with each element listing the
%         single unit IDs (same doubles as the ID's in the 'unitAll' property of EpyData) involved in the assembly. If no assembly is
%         found, the cell contains {NaN} for the concerned session.
%
%   Dependencies:
%     - compactSpikes: Removes gaps in single unit IDs for ICA compatibility.
%     - callISAC: Computes ISAC-based co-activation matrix (binary weights).
%
%   Requirements:
%     - The `unitsAll` property of the object must be pre-filled.
%     - `refDatesSessionTimestamps` must be populated with session start/end times.

%% Arguments block
arguments
    this (1,1) EpyData
    opt.rat (1,:) double = this.ratID
    opt.timeWindow2DetectAssemb (1,1) double = 0.03 % 30 ms time window (in seconds) for assembly detection
end

%% Compute

% Initialize storage
if isempty(this.assemblies)
    this.assemblies = cell(numel(this.ratID), 1);
end

% Loop on rats
for r = 1:numel(opt.rat)
    
    idxRat = find(this.ratID == opt.rat(r));
    spikes = this.unitsAll{idxRat,1};

    % Loop on the sessions of the rats
    for s = 1:numel(this.sessionPath{idxRat,1})

        % Initialisation
        this.assemblies{idxRat,1}.compositionISAC{s,1} = {NaN};

        % Get the time window specific to each session
        SessionTimeWindow = this.refDatesSessionTimestamps{idxRat,1}.cumulative(s,:);
        % Keep only spikes in the specified time window
        spikeTimes = spikes(:,1);
        selectedSpikes = spikes(spikeTimes > SessionTimeWindow(1) & spikeTimes < SessionTimeWindow(2), :);

        % Avoid having gap into single units ID whch will turn in "polluting false single units" for ISAC
        % Rename the single unit IDs for it to have no gap
        ICA_spikes = compactSpikes(selectedSpikes);
        singUnitIDTableRef = sort(unique(selectedSpikes(:,2)));

        % Compute assembly detection with ISAC
        tic;
        weights_ISAC = callISAC(ICA_spikes, opt.timeWindow2DetectAssemb); % dimension (neurons x assemblies)
        elapsedTime = toc;
        fprintf('ISAC detection executed in %.4f seconds.\n', elapsedTime);
        
        if any(weights_ISAC ~= 1 & weights_ISAC ~= 0)
            fprintf("You're a looser there is no only 0 and 1 in weights_ISAC")
        else  
            logical_ISAC = logical(weights_ISAC);
        end

        % Assemblies storage with single unit IDs which compose them
        if ~isempty(logical_ISAC)
            assemblies = cell(1,size(logical_ISAC,2));
            for a = 1:size(assemblies,2)
                assemblies{1,a} = singUnitIDTableRef(logical_ISAC(:,a));
            end
        else 
            assemblies = {NaN};
        end

        this.assemblies{idxRat,1}.compositionISAC{s,1} = assemblies;

    end

end