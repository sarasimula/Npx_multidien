function freqFiltered =  plotTestBestSignalProcessing(this, opt)
% plotTestBestSignalProcessing - Visualizes signal interpolation and filtering methods for IS/seizure frequency.
%
% Description:
%   ALL THE PARAMETERS INSIDE THE FUNCTIONS HAS BEEN TUNED AND ARE SPECIFICS TO TO opt.binCount = 3600s = 1s
%   This function retrieves the binned frequency of interictal spikes (IS) or seizures 
%   for a given rat, interpolates missing data using various methods (harmonics, spline, makima, pchip),
%   and filters the signal using a low-pass FIR filter to suppress fast fluctuations (e.g., below 12h).
%
% Inputs:
%   this : (1,1) EpyData instance.
%   opt  : structure with optional fields:
%       - rat      : IDs of the rats to include (default: all rats).
%       - nature   : "IS" or "seizures" (default: "IS").
%       - state    : "all", "sleep", or "wake" (default: "all").
%       - binCount : Time bin size in seconds (default: 3600 s = 1 hour).
%       - range    : Max time (s) to include (default: 0 = full range).
%
% Output:
%   freqFiltered : The filtered frequency signal (1D vector).
%
% Required functions : harmonicInterp

% Parse arguments
arguments
    this (1,1) EpyData
    opt.rat(1,:) double = this.ratID  % optionnal, by défault all the rats
    opt.nature (1,1) string {mustBeMember(opt.nature, ["IS", "seizures"])} = "IS"
    opt.state (1,1) string {mustBeMember(opt.state, ["all", "sleep", "wake"])} = "all" % By default, all IS regardless of the state
    opt.binCount (1,1) double = 3600  % Bin size in seconds for the count
    opt.range (1,1) double {mustBeNonnegative} = 0
end

for i = 1:numel(opt.rat)
    ratID = opt.rat(i);

    % Get event frequency over time
    freq = this.getFrequency("rat", ratID, "nature", "IS", ...
        "state", opt.state, "binCount", opt.binCount).freq{1}'; % * opt.binCount;
    t = (0:numel(freq)-1) / 24; % time in days

    % Assure to have raw lines
    freq(:);
    t(:);

    %% Smooth valid segments (without crossing NaN gaps) 
    % before doing the interpolation (for it to work better)
    windowSize = 15; % ajustable
    freqSmoothed = smoothSegments(freq, windowSize);
    
    %% Perform the interpolation (with different interpolations) only on the gap segments
    % Harmonics interpolation
    nHarmonics = 100;
    freqHarmonics = freq;
    freqInterpFull = harmonicInterp(t, freqSmoothed, nHarmonics);
    idxNaN = isnan(freq);
    freqHarmonics(idxNaN) = freqInterpFull(idxNaN); % At this step, freqHarmonics has filled gaps

    % Smooth the complete signal (there is no double smoothing) 
    freqHarmonicsSmooth = smoothdata(freqHarmonics, 'gaussian', 15); % The kernel must be < 18 (smallest segment between NaN blocks in the dataset)

    % Alternative interpolations
    % freqSpline = freq;
    % freqSpline(idxNaN) = spline(t(~idxNaN), freq(~idxNaN), t(idxNaN));
    freqSpline = spline(t, freqSmoothed, t);
    freqMakima = makima(t, freqSmoothed, t);
    freqPchip = pchip(t, freqSmoothed, t);

    %% Plot interpolated signals with offsets for readability
    offset = 0.1 * max(abs(freq), [], 'omitnan');  % 10% of max amplitude as offset
    
    figure; hold on;
    plot(t, freq, 'k.-');               % Plot real data

    % plot(t, freqSmoothed, 'g-');
    % plot(t, freqHarmonics, 'r-', 'LineWidth', 1.5);
    plot(t, freqHarmonicsSmooth, 'g-', 'LineWidth', 1.5)
    plot(t, freqSpline, 'b--', 'LineWidth', 1.5);
    plot(t, freqMakima, 'k--', 'LineWidth', 1.5);
    plot(t, freqPchip, 'r--', 'LineWidth', 1.5);
    xlabel('Temps [jours]'); ylabel('Amplitude ± offset');
    legend('Signal original', 'Interpolé (harmoniques)', 'Interpolé (spline)', 'Interpolé (makima)', "Interpolé (Pchip)");

    % FIR Filter design
    % Paramètres
    Fs = 24;                  % échantillons par jour
    Fc = 2;                   % fréquence de coupure (2 cycles/jour <=> 12h)
    Fp = 1.5;                     % bande passante (sous Fc)
    Fst = 2.5;                % bande stop (au-dessus de Fc)
    Ap = 1;                   % ondulation max (dB) dans bande passante
    Ast = 60;                 % atténuation min (dB) dans bande stop
    Fo = 60;                    % Filter order

    % Design du filtre passe-bas FIR
    d = designfilt('lowpassfir', ...
        'PassbandFrequency', Fp, ...
        'StopbandFrequency', Fst, ...
        'SampleRate', Fs, ...
        'DesignMethod', 'equiripple', ...
        'FilterOrder', Fo);
    % 'PassbandRipple', Ap, ...
    % 'StopbandAttenuation', Ast, ...


    % Apply zero-phase filtering
    freqFiltered = filtfilt(d, freqHarmonics);

    figure;
    hold on;
    plot(t, freq, 'k.-');
    plot(t, freqHarmonics + 5*offset, 'r-', 'LineWidth', 1.5);
    plot(t, freqFiltered, 'b-', 'LineWidth', 1.5);
    legend('Original', 'Interpolé (harmoniques)', 'Filtré (FIR, 12h)');
    xlabel('Temps [jours]'); ylabel('Amplitude ± offset');
    title(sprintf('Filtrage FIR 12h avec filtfilt - %s', ratID));
end
end


%        %        %         %        %        %         %        %        %         % 
function freqInterp = harmonicInterp(t, freqRaw, nHarmonics)
    % harmonicInterp - Fits missing values using a harmonic (Fourier-like) regression.
    %
    % Inputs:
    %   t          : Time vector (in days)
    %   freqRaw    : Input signal with NaNs
    %   nHarmonics : Number of harmonics in the regression (e.g. 5–100)
    %
    % Output:
    %   freqInterp : Interpolated signal (NaNs replaced)

    t = t(:);       
    freqRaw = freqRaw(:);

    % Keep only valid (non-NaN) data
    isValid = ~isnan(freqRaw);
    t_valid = t(isValid);
    y_valid = freqRaw(isValid);

    % Center time for numerical stability
    t0 = mean(t_valid);
    t_shift = t_valid - t0;
    t_all_shift = t - t0;

    % Design harmonic regression matrix
    X = ones(length(t_valid), 1);
    X_all = ones(length(t), 1);
    for k = 1:nHarmonics
        X = [X, cos(2*pi*1/4*k*t_shift), sin(2*pi*1/4*k*t_shift)];
        X_all = [X_all, cos(2*pi*1/4*k*t_all_shift), sin(2*pi*1/4*k*t_all_shift)];
    end

    % Least-squares fit
    beta = X \ y_valid;

    % Reconstruct full signal
    freqInterp = X_all * beta;
end

function freqSmoothed = smoothSegments(freqRaw, windowSize)
% smoothSegments - Smooths only valid continuous segments of a 1D signal separated by large enough NaN gaps.
%
% This function applies a Gaussian smoothing filter to continuous segments of a signal vector,
% where segments are defined as sequences of non-NaN values *delimited* by ≥5 consecutive NaNs.
% NaN regions and short segments (shorter than the smoothing window) are left untouched.
%
% Inputs:
%   freqRaw     : (Nx1 double) original signal vector that may contain NaNs
%   windowSize  : (scalar, optional) size of the smoothing window in samples (default: 10)
%
% Output:
%   freqSmoothed : (Nx1 double) signal with smoothed valid segments and untouched NaNs

    if nargin < 2
        windowSize = 10;  % Default smoothing window size
    end

    freqSmoothed = freqRaw;          % Initialize output
    isNaN = isnan(freqRaw);          % Logical index of NaNs
    minGap = 5;                      % Minimum number of consecutive NaNs to define a break

    % Identify NaN gaps of sufficient length (i.e., real segmentation points)
    isGap = false(size(freqRaw));   % Preallocate logical vector
    d = diff([false, isNaN, false]);  % Detect transitions into and out of NaN runs
    startNaN = find(d == 1);          % Indices where NaN sequences start
    endNaN = find(d == -1) - 1;       % Indices where NaN sequences end

    % Mark NaN sequences that are long enough to be considered actual gaps
    for k = 1:length(startNaN)
        if endNaN(k) - startNaN(k) + 1 >= minGap
            isGap(startNaN(k):endNaN(k)) = true;
        end
    end

    % Identify continuous segments between those gaps
    d = diff([true, isGap, true]);    % Detect transitions into and out of valid segments
    segStart = find(d == -1);         % Start indices of valid segments
    segEnd   = find(d == 1) - 1;      % End indices of valid segments

    % Smooth each valid segment individually
    for i = 1:length(segStart)
        idxStart = segStart(i);
        idxEnd   = segEnd(i);
        segment = freqRaw(idxStart:idxEnd);

        % Only smooth if the segment is long enough
        if numel(segment) >= windowSize
            freqSmoothed(idxStart:idxEnd) = smoothdata(segment, 'gaussian', windowSize);
        end
    end
end
