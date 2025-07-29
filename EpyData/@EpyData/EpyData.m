classdef EpyData
% EpyData Class to manage and structure multimodal electrophysiological data from chronic epilepsy models in rodents.
% 
% This class provides a unified interface for loading, organizing, and accessing key session-level data from two types 
% of acquisition systems: Neuropixels ("NPX") and telemetry ("TM"). It parses metadata and experimental events such as 
% interictal spikes (IS), seizures, sleep/wake intervals, and recording interruptions, aggregating them by animal (rat).
%
% Each row of the properties corresponds to a unique subject (rat), allowing for scalable batch processing across multiple
% sessions. The design ensures harmonized outputs independent of acquisition type, facilitating comparative or pooled analyses.
%
% Required functions : 
%   - FMAToolbox
%   - the helper functions : /mnt/hubel-data-103/Corentin/EpyData/Helper_Functions
%   - the Pietro's runBatch function : /mnt/cortex-data-319/Pietro/Code/MATLAB/Personal/runBatch.m
%
% Used files per rats :
%   - Neuropixel :
%       * .RealTimeCat.evt : 
%       * .cat.evt
%       * .cortex_IS_corr
%       (* .res)
%       (* .clu)
%       * .sws
%       * .seizures
%
%   - Telemetry :
%       * _sessions_total_points.json
%       * _N.sws
%       * _N.rem
%       * _N.IS
%       * .seizures
%
% Constructor :
%   obj = EpyData(PathBatch, DataType, opt)
%       - PathBatch : string. Path to the batch folder containing the session subfolders (I used a order batch by rats and days, it should work on an unordered batch but never tried).
%       - DataType  : string. Must be either "NPX" or "TM", to specify the acquisition system.
%       - opt       : optional struct with fields :
%           * ZT             : double (default = NaN). Zeitgeber Time reference (used for NPX only).
%           * tsMicroAwake   : double (default = 0). Threshold below which micro-awakenings intervals are merged into sleep intervals.
%           * verbose        : boolean (default = true). Display processing steps in command window.
%           * units          : boolean (default = false). Also add units spikes timestamps in epydata.
%
% Properties:
%   sessionPath            - Cell array of strings : each row contains paths to all folders (sessions) associated with one rat.
%   ratID                  - Numeric array: unique identifier for each rat (ex: 596, 597, 1214...).
%   recDates               - Nx2 datetime: start and end of recording in real-world calendar time.
%   refDatesTimestamps     - Nx2 datetime: time reference of all events, possibly shifted to Zeitgeber Time (ZT).
%   
%   sleepInt               - Cell array of Nx2 doubles: timestamps (s) of sleep intervals (SWS + REM).
%   swsInt                 - Cell array of Nx2 doubles: timestamps (s) of slow waves sleep intervals (SWS).
%   unrecInt               - Cell array of Nx2 doubles: timestamps (s) of recording interruptions.
%   wakeInt                - Cell array of Nx2 doubles: awake intervals deduced from sleep and unrecorded intervals.
%   ISTimestamps           - Cell array of Nx1 vectors: timestamps (s) of interictal spike events.
%   seizuresTimestamps     - Cell array of Nx2 doubles: start and end timestamps (s) of seizure events.
%   unitsAll               - Cell array of Nx2 doubles: concerns ONLY good single unit & computed ONLY with opt. 
%                               First column : timestamp (in seconds) of spikes
%                               Second column : cumulative (across sessions) single unit ID.
%   ZT                     - Double between 0 and 24 or NaN : Zeitgeber Time offset used for temporal alignment (for Neuropixel). 
%                          Normally ZT = 7 (7 am) corresponds to light on. NaN if no rescaling to ZT.
%   samplingFreq           - Sampling frequency (Hz) of the signal (30kHz for NPX, 512Hz for TM).
%   rangeCutIS             - Numeric array : maximum timestamp in seconds until correct and exploitable recording. Created thanks to computeRangeCutIS.
%   PSDreal & PSDsurro     - Cell array of structure having same fileds for all states (.all, .sleep, .wake). Empty at the contruction. Created thanks to computePSD. 
%                               Each state field is composed of following . 
%                                   - .f : vector of grid of frequencies used for the PSD computation (different for TM and NPX). Should now be the same grid for the real data & PSD associated.
%                                   - .pxx : vector of not normalized power values of the PSD on real IS count serie (resp. (nSurrogates x 1) cell array of PSD of shuffles).
%   multidianCycl          - Empty at the construction. Created thanks to analyzePSD methods.
%                            Cell array (each rat) of structure having fields :
%                               - f0 : Significant peak frequencies [Hz] in the multidian range (3 to 15 days, or up to the duration of the session) where real PSD is significantly higher than surrogates (computed using the 'all' state).
%                                  Ordered in ascendant order of importance of the peak. 
%                               - importanceFact_f : Relative importance of the peak associated to the frequency of f0.
%                                  Ordered in the same order as f0.
%                               - fWidth : Boundaries frequency [fStart, fEnd] describing each peak at f0. Calculated finding the f associated to power_f0/2.
%                               - sigBandEdges :  Frequency bins [fStart, fEnd] significantly above surrogates.
%                               - freqISFiltered : structure with field of all state (.all, .sleep, .wake) giving the multidien temporal serie rhythm (after band passing at each f0).
%                               - incrFromShufSleep : Difference between normalized real PSD of IS count and surrogate mean at the different f0, in the 'sleep' state
%                               - incrFromShufWake : Same as above, for the 'wake' state
%                               - areaFromShufSleep : At every f0 peak, area under the normalized real PSD and the nomalized surrogates in the range of the peak width, in the 'sleep' state
%                               - areaFromShufWake : Same as above, for the 'wake' state
%   assemblies             - Empty at the construction. Created thanks to getCellAssemblies and completed by activationCellAssemblies. 
%                            Cell array (each rat) of structure having fileds :
%                               - compositionISAC: Each cell represents a session and contains a (T x 2) matrix where T is the time activation 
%                                 of the assembly which ID is given in the associated line in column 2.  
%                               - activationISAC : Each cell represents a session and contains a 1xA cell array (A = number of assemblies), with each element listing the
%                                 single unit IDs (same doubles as the ID's in the 'unitAll' property of EpyData) involved in the assembly. If no assembly found, cell contains {NaN}.
%
% Methods:
%   COMPUTATION :
%       IS :                    - Retrieves timestamps of interictal spikes (IS) for selected rats and brain states.
%       stateFractionPerBins :  - Computes the time spent in a specific state per time bin for one rat.
%       getFrequency :          - Computes time-binned frequency of events (IS or seizures) per rat.
%       surrogateIS :           - Generates surrogate interictal spike (IS) frequency time series.
%       computePSD :            - Computes Lomb-Scargle power spectral density (PSD) of interictal spike (IS) count data
%       analyzePSD :            - Analyzes Power Spectral Density and identify significant multidian bands.
%       getStateIntervals :     - Returns time intervals corresponding to a given brain state for selected rats.
%       ect ... (may check the log for more informations)
%
%   PLOT :
%       plotComparativeViolin : - Compare two group distributions with violin plot, box plots and scatters (+ stats).
%       plotPSDPerRat         : - Plots IS time series and PSDs (real + surrogates) per state. A figure per rat.
%       ect ... (may check the log for more informations)
%
% Author : Corentin 2025 @Copyright
  

properties (GetAccess = public, SetAccess = protected)
  sessionPath (:,1) cell
  ratID (:,1) double
  recDates (:,2) datetime
  refDatesTimestamps (:, 2) datetime
  refDatesSessionTimestamps (:,1) cell % each element is a structure
  sleepInt (:,1) cell
  swsInt (:,1) cell
  unrecInt (:,1) cell
  wakeInt (:,1) cell
  ISTimestamps (:,1) cell
  seizuresTimestamps (:,1) cell
  ZT (1,1) 
  unitsAll (:,1) cell
  samplingFreq (1,1) double
  PSDreal (:,1) cell        % Cell of structure with all the states for real PSD
  PSDsurro (:,1) cell       % Cell of structure with all the states for surrogates PSD
  multidienCycl (:,1) cell  % Each element is a structure containing the three state "all", "sleep" and "wake". For each state a f0 is saved if a significative frequency for multidien rhythm is found. The amplitude (A) and phase at t=0 (phi_t0) of the best sinusoïde is also kept here. To reconstruct
  rangeCutIS (:,1) double   % For IS, timestamps in seconds (from the refDatesTimestamps beginning) where the signal should be cut not to have edge effect (commonly used by default as opt.range)
  assemblies (:,1) cell     % Cell of structure with for each session the cell assemblie composition using PCA method (different weight) and ISAC method (no weight).
end
  
  methods  (Access = public)

    function obj = EpyData(PathBatch,DataType, opt)
        %       - PathBatch : string. Path to the batch folder containing the session subfolders (I used a order batch by rats and days, it should work on an unordered batch but never tried).
        %       - DataType  : string. Must be either "NPX" or "TM", to specify the acquisition system.
        %       - opt       : optional struct with fields :
        %           * ZT             : double (default = NaN). Zeitgeber Time reference (used for NPX only).
        %           * tsMicroAwake   : double (default = 0). Threshold below which micro-awakenings intervals are merged into sleep intervals.
        %           * verbose        : boolean (default = true). Display processing steps in command window.
        %           * units          : boolean (default = false). Also add units spikes timestamps in epydata

      arguments
          PathBatch (1,1) string
          DataType (1,1) string {mustBeMember(DataType, ["NPX", "TM"])}
          opt.ZT (1,1) double = NaN
          opt.tsMicroAwake (1,1) double = 0
          opt.verbose (1,1) {mustBeMember(opt.verbose, [0,1])} = true
          opt.units (1,1) {mustBeMember(opt.units, [0,1])} = false
      end
      
      % Two different pipelines depending on the acquisition system.
      switch DataType
          
          case "NPX"
              obj = obj.LoadNPX(PathBatch, opt);
          case "TM"
              obj = obj.LoadTM(PathBatch, opt);
      end

      % Compute the cutoff range timestamps for the instances
      obj = obj.computeRangeCutIS;

    end

  end

    methods (Access = private)
        
        function this = LoadNPX(this, PathBatch, opt)
            % LoadNPX Load and structure data from Neuropixels sessions.
            %
            % This method processes a batch of Neuropixels sessions by calling specific
            % helper functions (e.g., getMetadataNPX, getEpilepticEventsNPX), aggregating
            % data by rat. It applies optional Zeitgeber Time correction and sleep merging.
            %
            % Inputs:
            %   - PathBatch : string. Path to the folder containing Neuropixels sessions.
            %   - opt       : struct. Optional parameters:
            %       * ZT             : double. Zeitgeber Time for temporal alignment (default = NaN).
            %       * tsMicroAwake   : double. Threshold (s) for merging micro-awakenings (default = 0).
            %       * verbose        : boolean. Display processing logs (default = true).
            %       * units          : boolean. Get session units (default
            %
            % Output:
            %   - Updated object with EpyData properties for Neuropixels data.

            arguments
                this EpyData
                PathBatch (1,1) string
                opt struct 
            end

            if opt.verbose 
                fprintf("Getting data for each day session \n") 
            end
            [sessionPath, ratID, date, realTime] = runBatch(PathBatch, @getMetadataNPX, verbose = opt.verbose);
            args = {opt.tsMicroAwake, opt.ZT, opt.verbose, opt.units};
            if opt.units
                [IS_rescale, seizures_rescale, sleep_rescale, sws_rescale, wake_rescale, ~, ~, units_rescale] = runBatch(PathBatch, @getEpilepticEventsNPX, args, verbose = opt.verbose);
            else
                [IS_rescale, seizures_rescale, sleep_rescale, sws_rescale, wake_rescale, ~, ~, ~] = runBatch(PathBatch, @getEpilepticEventsNPX, args, verbose = opt.verbose);
                units_rescale = {};
            end

            if opt.verbose 
                fprintf("Aggregating and resacling data from sessions by rat \n")
            end
            [this.ratID, this.recDates, this.sessionPath] = aggregateMetadataNPX(sessionPath, ratID, date, realTime);
            [this.ISTimestamps, this.seizuresTimestamps, this.sleepInt, this.swsInt, this.unrecInt, this.wakeInt, this.refDatesSessionTimestamps, ...
                this.refDatesTimestamps, this.ZT, this.unitsAll] = aggregateEpylepticDataNPX(ratID, IS_rescale, seizures_rescale, sleep_rescale, sws_rescale, wake_rescale, units_rescale); % units_rescale is optional

            % Update the last element of the properties
            this.recDates(:,2) = this.refDatesTimestamps(:,2); % give the true end value of RefDatesTimestamp to RecDates
            this.samplingFreq = 30000;
        end

        function this = LoadTM(this,PathBatch, opt)
            % LoadTM Load and structure data from Telemetry sessions.
            %
            % This method processes a batch of Telemetry sessions by calling specific
            % helper functions (e.g., getSessionMetadataTM, getStateTM, getIsTM), aggregating
            % data by rat. No Zeitgeber Time alignment or recording interruption detection
            % is computed for TM data and sleep merging is optional.
            %
            % Inputs:
            %   - PathBatch : string. Path to the folder containing Telemetry sessions.
            %   - opt       : struct. Optional parameters:
            %       * tsMicroAwake   : double. Threshold (s) for merging micro-awakenings (default = 0).
            %       * verbose        : boolean. Display processing logs (default = true).
            %
            % Output:
            %   - Updated object with EpyData properties for Telemetry data.
            arguments
                this EpyData
                PathBatch (1,1) string
                opt struct
            end

            this.samplingFreq = 512; % EEG frequency

            if opt.verbose
                fprintf("Getting MetaData \n")
            end
            [this.sessionPath, ratID, recDates, cumulativeTimes] = runBatch(PathBatch, @getSessionMetadataTM, verbose = opt.verbose);
            
            this.ratID = cell2mat(ratID);
            this.recDates = vertcat(recDates{:});
            refDatesSessionTimestamps = cellfun(@(x) [x(1:end-1), x(2:end)], cumulativeTimes, 'UniformOutput', false);
            for i=1:numel(ratID)
                this.refDatesSessionTimestamps{i,1}.cumulative = refDatesSessionTimestamps{i,1};
                this.refDatesSessionTimestamps{i,1}.date = repmat(this.recDates(i,1),size(refDatesSessionTimestamps{i,1},1), size(refDatesSessionTimestamps{i,1},2)) + seconds(refDatesSessionTimestamps{i,1});
                
                % Sanity Check
                if this.refDatesSessionTimestamps{i,1}.date(end,2) ~= this.recDates(i,2)
                    error("Try again")
                end
            end
            
            % No Zeit Geber disponible for telemetry
            this.ZT = NaN; 
            this.refDatesTimestamps =  this.recDates;%
            
            % No interruption in the recording
            this.unrecInt = cell(numel(this.ratID),1);

            if opt.verbose
                fprintf("Getting state intervals \n")
            end
            args = {opt.tsMicroAwake, opt.verbose};
            [this.sleepInt, this.swsInt, this.wakeInt] = runBatch(PathBatch, @getStateTM, args, verbose = opt.verbose);

            if opt.verbose
                fprintf("Getting interictal spikes timestamps \n")
            end
            args = {opt.verbose};
            [this.ISTimestamps] = runBatch(PathBatch, @getIsTM, args, verbose = opt.verbose);

            if opt.verbose
                fprintf("Getting seizures timestamps \n")
            end
            args = {opt.verbose};
            this.seizuresTimestamps = runBatch(PathBatch, @loadSeizuresTM, args, verbose = opt.verbose);

        end
    end


end