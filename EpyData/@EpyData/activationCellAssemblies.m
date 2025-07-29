function this = activationCellAssemblies(this, opt)
%ACTIVATIONCELLASSEMBLIES Computes temporal activation of detected cell assemblies
%
%   This method computes the temporal activation of previously detected cell assemblies
%   (stored in `compositionISAC`) for each rat and session. It uses the ISAC approach
%   on binned spike trains and projects activity onto the assembly templates to quantify
%   their activation over time. Activations are stored in the `activationISAC` field of the
%   `assemblies` property.
%
%   Syntax:
%     this = activationCellAssemblies(this, opt)
%
%   Inputs:
%     this                             - (1x1 EpyData) Instance of EpyData class, must contain
%                                        a non-empty `unitsAll` property and precomputed assemblies
%                                        in `compositionISAC`.
%
%     opt (Name-Value Arguments):
%       .rat                           - (1xN double) IDs of rats to process (default: all rats in `this.ratID`)
%       .timeWindow2ActivAssemb        - (1x1 double) Size of the sliding time window [in seconds] 
%                                        for computing activation (default: 0.03 s = 30 ms)
%
%   Changed Properties:
%     this.assemblies{r}.activationISAC{s} 
%       - (double matrix) Each cell contains a (T x 2) matrix where T is the time activation 
%         of the assembly which ID is given in the associated line in column 2. 
%
%   Dependencies:
%     - compactSpikes: Remaps spike unit IDs to continuous range for matrix operations.
%     - ICActivations: Computes activation time series for assemblies from spike trains and 
%                      assembly templates.
%
%   Requirements:
%     - The 'unitsAll' property must be populated.
%     - 'compositionISAC' must be precomputed via `getCellAssemblies`.

%% Arguments block
arguments
    this (1,1) EpyData
    opt.rat (1,:) double = this.ratID
    opt.timeWindow2ActivAssemb (1,1) double = 0.03 % 30 ms time window (in seconds) for assembly activation
end
%%

% Loop on rats
for r = 1:numel(opt.rat)

    idxRat = find(this.ratID == opt.rat(r));
    spikes = this.unitsAll{idxRat,1};  
    % ratMatrices = {};  % init

    % Loop on the sessions of the rats
    for s = 1:numel(this.sessionPath{idxRat,1}) 

        % Initialisation
        this.assemblies{idxRat,1}.activationISAC{s,1} = NaN;

        compo = this.assemblies{idxRat,1}.compositionISAC{s,1};

        isNaNCell = iscell(compo) && isequal(size(compo), [1 1]) ...
             && isnumeric(compo{1}) && isscalar(compo{1}) && isnan(compo{1});


        if ~isNaNCell
            %% Rename spikes to be continuous ids
            % Get the time window specific to each session
            SessionTimeWindow = this.refDatesSessionTimestamps{idxRat,1}.cumulative(s,:);
            % Keep only spikes in the specified time window
            spikeTimes = spikes(:,1);
            selectedSpikes = spikes(spikeTimes > SessionTimeWindow(1) & spikeTimes < SessionTimeWindow(2), :);

            % Avoid having gap into single units ID whch will turn in "polluting false single units" for ICA
            % Rename the single unit IDs for it to have no gap
            ICA_spikes = compactSpikes(selectedSpikes);
            nSingUnit = max(unique(ICA_spikes(:,2)));
            singUnitIDTableRef = sort(unique(selectedSpikes(:,2)));


            %% Reconstruct the matrix output of callISAC function
            sessionCellCompAssemblies = this.assemblies{idxRat}.compositionISAC{s};
            % Convert each stored assembly (list of unit indices) to binary vector
            sessionMatrixCompAssemblies = cell2mat(cellfun(@(x) reshape(accumarray( arrayfun(@(y) find(singUnitIDTableRef == y), x),1,[nSingUnit 1]),[],1), ...
                    sessionCellCompAssemblies, 'UniformOutput', false));

            %% Compute activation
            fprintf('ISAC activation started for session %d of rat %d. \n', s, opt.rat(r));
            tic;
            try
                ISAC_activations = ICActivations(ICA_spikes, sessionMatrixCompAssemblies, opt.timeWindow2ActivAssemb);

                % Catch required only for session 2 and session 5 of rat 794
            catch ME
                % Check if the error is due to memory allocation
                if contains(ME.message, 'array exceeds maximum array size preference', 'IgnoreCase', true)
                    warning('Memory allocation error for session %d of rat %d : ICActivations on several fragments of the spikes', s, opt.rat(r));
                    % What raises an error
                    maxElemAllowed = floor(503.4e9 / 8); % 8 bits per element and maximum size =  503.4 GB
                    
                    % Number of intervals to used ICActivations
                    nIntervals = size(Bins(round(ICA_spikes(1,1)),ICA_spikes(end,1),opt.timeWindow2ActivAssemb,opt.timeWindow2ActivAssemb/2),1);

                    nElements = nIntervals*size(sessionMatrixCompAssemblies,2); 

                    % Minimum number of time segments that should be used
                    % for ICActivation to not raise the allocation error
                    nSegments = ceil(nElements/maxElemAllowed);

                    edges = linspace(ICA_spikes(1,1),ICA_spikes(end,1), nSegments+1);
                    ISAC_activations = [];
                    for seg = 1:nSegments
                        ICA_spikes_seg = ICA_spikes(ICA_spikes(:,1) > edges(seg) & ICA_spikes(:,1) < edges(seg+1), :);
                        ISAC_activations_seg = ICActivations(ICA_spikes_seg, sessionMatrixCompAssemblies, opt.timeWindow2ActivAssemb);
                        ISAC_activations = [ISAC_activations; ISAC_activations_seg];
                    end
    
                else
                    rethrow(ME);  % Propagate other unexpected errors
                end
            end
            elapsedTime = toc;
            fprintf('ISAC activation executed in %.4f seconds.\n', elapsedTime);

            %% Store
            this.assemblies{idxRat,1}.activationISAC{s,1} = ISAC_activations;
        end

        % ratMatrices{end+1} = sessionMatrix;  % accumulate
    end

    % Combine all into block-diagonal matrix
    % ratMatrixCompAssemblies = blkdiag(ratMatrices{:});

    

end
