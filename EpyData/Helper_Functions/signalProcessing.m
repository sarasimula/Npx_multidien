function [freqFiltered, freqs, A] =  signalProcessing(freqRaw, opt)
% signalProcessing - Interpolates and filters a frequency signal using harmonic regression and FIR filtering.
%
% syntax :
%   [freqFiltered, freqs, A] =  signalProcessing(freqRaw, nHarmonics, opt)
%
% Description:
%   THE PARAMETERS BEFORE THE FILTER HAS BEEN TUNED AND ARE SPECIFICS TO opt.binCount = 3600s = 1H
%   This function takes a raw frequency signal (e.g., binned interictal spike),
%   interpolates missing values using harmonic regression with `nHarmonics` sine and cosine terms, 
%   and applies a low-pass equiripple FIR filter to remove fast components (by default under 12 h cycles).
%
% Inputs:
%   freqRaw      : (1,N) vector of raw frequency values (NaNs allowed for missing data).
%   nHarmonics   : (scalar) number of harmonics used for the interpolation (default: 100).
%   opt          : structure with optional fields for filtering:
%       - Fs   : sampling rate (samples per day), default = 24 (i.e., 1 point per hour).
%       - Fc   : cutoff frequency in cycles/day (informative, not directly used here).
%       - Fp   : passband frequency in cycles/day (e.g., 1.5).
%       - Fst  : stopband frequency in cycles/day (e.g., 2.5).
%       - Fo   : filter order (e.g., 60).
%
% Outputs:
%   freqFiltered - interpolated and filtered version of freqRaw
%   freqs        - frequencies of each harmonic (cycles/day)
%   A            - amplitude of each harmonic component (proxy for power)
%
% Dependencies:
%   - Requires Signal Processing Toolbox (for designfilt and filtfilt).
%
% Uses the help function : smoothSegments

% Parse arguments
arguments
    freqRaw (1,:)
    opt.nHarmonics (1,1) double {mustBeInteger} = 100 % Not a problem if it overfits the high frequency noise because lowfilter pass after
    opt.Fs (1,1) double = 24;        % samples per day
    opt.Fc (1,1) double = 2;         % cutoff frequency (2 cycles/day <=> 12h period)
    opt.Fp (1,1) double = 1.5;       % passband frequency (below Fc)
    opt.Fst (1,1) double = 2.5;      % stopband frequency (above Fc)
    opt.Fo (1,1) double = 60;        % filter order
end

% Time vector in days, assuming uniform hourly bins (Fs = 24 samples/day)
t = (0:numel(freqRaw)-1) / opt.Fs; 
t = t(:);                           % Ensure column vector
freqRaw = freqRaw(:);              % Ensure column vector

% Smooth valid segments (without crossing NaN gaps) before doing the interpolation (for it to work better)
windowSize = 15; % ajustable
freqSmoothed = smoothSegments(freqRaw, windowSize);

% Perform the interpolation (harmonic regression) to fill the gap
nHarmonics = 100;
freqHarmonics = freqRaw;
[freqInterpFull, beta] = harmonicInterp(t, freqSmoothed, nHarmonics);
idxNaN = isnan(freqRaw);
idxZeros = freqRaw==0;
freqHarmonics(idxNaN) = freqInterpFull(idxNaN); freqHarmonics(idxZeros) = freqInterpFull(idxZeros); % At this step, freqHarmonics contains freqRaw and has filled gaps (of NaN and zeros)

% Design low-pass FIR filter with equiripple method
d = designfilt('lowpassfir', ...
    'PassbandFrequency', opt.Fp, ...
    'StopbandFrequency', opt.Fst, ...
    'SampleRate', opt.Fs, ...
    'DesignMethod', 'equiripple', ...
    'FilterOrder', opt.Fo);

% Apply zero-phase digital filtering (no phase distortion)
freqFiltered = filtfilt(d, freqHarmonics);

% Compute amplitudes specific to every frequency of the sinusoids
nHarm = (length(beta) - 1) / 2;
A = zeros(nHarm, 1);
for k = 1:nHarm
    idx_cos = 2*k;
    idx_sin = 2*k + 1;
    A(k) = sqrt(beta(idx_cos)^2 + beta(idx_sin)^2);
end
freqs = (1:nHarm)' * (2/5); % fréquence en cycles/jour, the 2/5 factor is coming from the 2/5 factor used for the harmonic matrix construction
% 2/5

end


%%        %        %         %        %    HELP FUNCTIONS    %         %        %        %         % 

function [freqInterp, beta] = harmonicInterp(t, freqRaw, nHarmonics)
    % harmonicInterp - Fits missing values using a harmonic (Fourier-like) regression.
    %
    % Inputs:
    %   t          : Time vector (in days)
    %   freqRaw    : Input signal with NaNs
    %   nHarmonics : Number of harmonics in the regression (e.g. 5–100)
    %
    % Output:
    %   freqInterp : Interpolated signal (NaNs replaced)
    %   beta       : the vector of regression coefficients for the harmonic (Fourier-like) model used to interpolate the signal.

    t = t(:);       
    freqRaw = freqRaw(:);

    % Identify valid (non-NaN) data points
    isValid = ~isnan(freqRaw);
    t_valid = t(isValid);
    y_valid = freqRaw(isValid);

    % Center time around mean for numerical stability
    t0 = mean(t_valid);
    t_shift = t_valid - t0;
    t_all_shift = t - t0;

    % Construct harmonic basis matrices (cosine and sine terms up to nHarmonics)
    X = ones(length(t_valid), 1);
    X_all = ones(length(t), 1);
    for k = 1:nHarmonics
        X = [X, cos(2*pi*2/5*k*t_shift), sin(2*pi*2/5*k*t_shift)];
        X_all = [X_all, cos(2*pi*2/5*k*t_all_shift), sin(2*pi*2/5*k*t_all_shift)];
    end

    % Solve least-squares regression: estimate harmonic coefficients
    warningState = warning('off', 'MATLAB:rankDeficientMatrix');
    beta = X \ y_valid;
    warning(warningState);

    % Reconstruct full interpolated signal
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
    d = diff([false, isNaN', false]);  % Detect transitions into and out of NaN runs
    startNaN = find(d == 1);          % Indices where NaN sequences start
    endNaN = find(d == -1) - 1;       % Indices where NaN sequences end

    % Mark NaN sequences that are long enough to be considered actual gaps
    for k = 1:length(startNaN)
        if endNaN(k) - startNaN(k) + 1 >= minGap
            isGap(startNaN(k):endNaN(k)) = true;
        end
    end

    % Identify continuous segments between those gaps
    d = diff([true, isGap', true]);    % Detect transitions into and out of valid segments
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
