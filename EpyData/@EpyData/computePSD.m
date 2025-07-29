function [this, PSDReal] = computePSD(this, opt)
% COMPUTEPSD - Compute Lomb-Scargle power spectral density (PSD) of interictal spike (IS) count data
%              using MATLAB's built-in `plomb` function. Computes PSD for real IS counts (returned)
%              and surrogate series (stored in the object because of a long time of computation).
%
% Syntax:
%   [this, PSDReal] = this.computePSD(opt)
%
% Description:
%   Computes the PSD of interictal spike (IS) counts using the Lomb-Scargle method, allowing for
%   uneven sampling or missing values. The real PSD is stored, and surrogate PSDs are optionally
%   generated and stored. Surrogates are recomputed only if not present or if 'force' is set to true.
%
%   Frequency resolution is adapted to the data type (Neuropixels or Telemetry) and uses a
%   logarithmic frequency grid.
%
% Modifications to `this` (EpyData object):
%   - Adds/updates the field `this.PSDreal{idx}.(state)` with:
%       • .f   : frequency vector [Hz]
%       • .pxx : PSD of the real IS count signal
%   - Adds/updates the field `this.PSDsurro{idx}.(state)` with:
%       • .f   : frequency vector [Hz] (shared with real)
%       • .pxx : cell array of PSDs for surrogate IS signals
%
% Inputs:
%   this        - (EpyData) Object instance containing IS data.
%   opt         - (struct, optional) Name-value options:
%       • rat        : [1×N double] Rat IDs to include (default = all rats in object)
%       • state      : ["all"|"sleep"|"wake"] State to filter IS (default = "all")
%       • binCount   : (double) Bin size in seconds (default = 3600)
%       • range      : (double) Duration of data to keep in seconds, [] = all (default = [])
%       • nSurro     : (int) Number of surrogate series to generate (default = 100)
%       • binShuffle : (double) Bin size for surrogate shuffling in seconds (default = 86400)
%       • force      : (logical) Force recomputation of surrogates even if results exist (default = false)
%
% Outputs:
%   this        - Updated EpyData object with PSD fields.
%   PSDReal     - Cell array [nRats×1] of structs with:
%                   • .f   : Frequency grid [Hz]
%                   • .pxx : PSD of real IS count signal
%
% Dependencies:
%   - surrogatesIS (method of EpyData)


arguments
    this (1,1) EpyData
    opt.rat(1,:) double = this.ratID
    opt.state (1,1) string {mustBeMember(opt.state, ["all", "sleep", "wake"])} = "all"
    opt.binCount (1,1) double = 3600                        % Time bin size for IS count in seconds
    opt.range (:, 1) double {mustBeNonnegative} = []        % Duration to keep, [] = all
    opt.nSurro (1,1) double {mustBeNonnegative, mustBeInteger} = 100
    opt.binShuffle (1,1) double {mustBePositive} = 86400    % Bin size for shuffling (surrogates) in seconds
    opt.force (1,1) logical = false
end

nRats = numel(opt.rat);                                   % Number of rats to process
PSDReal = cell(nRats, 1);                                 % Output: PSD of real data for each rat

% Define common frequency grid (log scale)
fmax = 1 / opt.binCount;                                  % Max frequency based on bin size
fminCommun = 1 / (3600*24*15);                            % Min frequency shared for our both population (period = 15 days)
NfCommun = 5000;                                          % Number of frequency points of the PSD on the shared frequency grid
FGridCommun = logspace(log10(fminCommun), ...             % Frequency grid (log scale) for the PSD
    log10(fmax), NfCommun);  

% Adapt frequency range depending on acquisition system
switch this.samplingFreq
    case 512                % TM
        fminTM = 1 / (3600*24*21);                                  % Min freq for TM (21-day period)
        NfTM =  floor( NfCommun * abs(log10(fminCommun)-log10(fminTM))/ ...
            abs(log10(fmax)-log10(fminCommun)));                    % Number of points for the extra grid
        FGrid = [logspace(log10(fminTM), log10(fminCommun), NfTM) , ...
            FGridCommun(2:end)];                                    % Concatenation of frequency grid (log scale) for the PSD
    case 30000              % NPX
        % Update fminCommun if fminNPX needs to be changed. 
        % In any case, we absolutely need fminTM <= fminNPX.
        FGrid = FGridCommun;
end

% Initialize storage for real and surrogates if not yet present
if isempty(this.PSDsurro)
    this.PSDsurro = cell(numel(this.ratID), 1);
end
if isempty(this.PSDreal)
    this.PSDreal = cell(numel(this.ratID), 1);
end

% Get real IS counts
ISReal = this.surrogatesIS("rat", opt.rat, "state", opt.state, ...
    "binCount", opt.binCount, "range", opt.range, "nSurro", 0);

% Get surrogates if needed
if opt.nSurro ~= 0
    ISSurro = this.surrogatesIS("rat", opt.rat, "state", opt.state, ...
        "binCount", opt.binCount, "range", opt.range, ...
        "nSurro", opt.nSurro, "binShuffle", opt.binShuffle);
end

% Loop through each rat
for j = 1:nRats
    ratIDj = opt.rat(j);
    ratIdx = find(this.ratID == ratIDj, 1); % Find rat index in object

    if isempty(ratIdx)
        warning("Rat ID %d not found in class. Skipping.", ratIDj);
        continue
    end

    % Create time vector corresponding to IS counts
    nBinsReal = numel(ISReal{j});
    t = this.refDatesTimestamps(ratIdx,1) + seconds(0:opt.binCount:(nBinsReal-1)*opt.binCount);

    % Compute PSD for real IS series
    [psdR, fR] = plomb(ISReal{j}(:), t, FGrid, 'psd');
    PSDReal{j} = struct("pxx", psdR, "f", fR); % Store result
    this.PSDreal{ratIdx}.(opt.state) = struct("pxx", psdR, "f", fR); 

    % Ensure structure exists for that rat
    if isempty(this.PSDsurro{ratIdx})
        this.PSDsurro{ratIdx} = struct();
    end

    % Check existing surrogate computation
    hasSurro = isfield(this.PSDsurro{ratIdx}, opt.state);
    if hasSurro && ~opt.force
        fprintf("Some PSD surrogates were already computed for rat %d (%s). Use of those former surrogates.\n To force the computation use the input force = true. \n", ...
            ratIDj, opt.state);
    end

    % Force computation or no previous result: (re)compute surrogate PSD
    if (~hasSurro && opt.nSurro ~= 0) || opt.force 
        fprintf("Computation of PSD surrogates for rat %d (%s). \n", ratIDj, opt.state);

        
        this.PSDsurro{ratIdx}.(opt.state) = struct();
        this.PSDsurro{ratIdx}.(opt.state).f = FGrid';

        nSurro = size(ISSurro{j}, 2); % Number of surrogates
        nBinsSurro = size(ISSurro{j}, 1);
        
        for k = 1:nSurro
            % Create time vector corresponding to IS counts
            tSurro = this.refDatesTimestamps(ratIdx,1) + seconds(0:opt.binCount:(nBinsSurro-1)*opt.binCount);
            % Lomb-Scargle periodogram of the zscored signal
            [this.PSDsurro{ratIdx}.(opt.state).pxx{k,1}] = plomb(ISSurro{j}(:,k), tSurro, FGrid, 'psd');
        end
    end
end
end