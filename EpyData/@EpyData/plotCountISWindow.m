function plotCountISWindow(this, opt)
% plotCountISWindow - Zoom the IS time series over a time window & 
%                   superpose Light On/Light Off thanks to ZeitGebber.
%
% syntax :
%   this = plotCountISWindow(this, opt)
%
% Inputs:
%   this : (1,1) EpyData object containing PSD and event data
%   opt  : struct with optional parameters:
%       - rat                : vector of rat IDs to process (default: all rats)
%       - state              : string - 'all' | 'sleep' | 'wake' (default = 'all')
%       - binCount           : bin size (s) for temporal binning (default: 3600)
%       - rangeWindow        : beggining and ending timestamp (s), for the zoom window
%
% Output:
%   this : updated EpyData object with potentially recomputed PSDs
%
% Methods and functions used : 
%   - getFrequency
%   - signalProcessing

arguments
    this (1,1) EpyData
    opt.rat(1,1) double = this.ratID(1)
    opt.state (1,1) string {mustBeMember(opt.state, ["all", "sleep", "wake"])} = "all" % By default, all IS regardless of the state
    opt.binCount (1,1) double = 3600
    opt.rangeWindow (:, 2) double {mustBeNonnegative} = [0, seconds(this.refDatesTimestamps(1,2) - this.refDatesTimestamps(1,1))]
end

idx = find(this.ratID == opt.rat);

freq = this.getFrequency("rat", opt.rat, "nature", "IS", ...
    "state", opt.state, "binCount", opt.binCount).freq{1}'; % * opt.binCount;
freqFiltered = signalProcessing(freq, nHarmonics=20);
t = (0:numel(freq)-1) / 24; % time in days


figure('Name', sprintf('Rat %d - Zoomed IS count', opt.rat), ...
    'Units', 'pixels', 'Position', [0, 0, 2400, 1600], ...
    'Visible','on');
hold on; 

% Plot real frequency
plot(t, freq, 'Color', [0 0 1], 'LineWidth', 1.5, 'DisplayName', 'IS count');

% Plot multidien rythm or just the smooth data series (in case of non multidien rhythm)
if isfield(this.multidienCycl{idx,1}, "freqISFiltered") && ~isempty(this.multidienCycl{idx,1}.freqISFiltered)
    plot(t, this.multidienCycl{idx,1}.freqISFiltered.(opt.state), 'Color', [0.1 0.2 0.1], 'LineWidth', 2.5, 'DisplayName', 'Multidien Rythm')
else
    plot(t, freqFiltered, 'Color', [1 0.5 0], 'LineWidth', 2.5, 'DisplayName', 'Smooth ')
end

% Seizures
if ~isempty(this.seizuresTimestamps{idx})
    seiz = this.seizuresTimestamps{idx};
    tEnd = seconds(this.refDatesTimestamps(idx,2) - this.refDatesTimestamps(idx,1));

    overlaySeizures(seiz, tEnd, numberHours = 6, cumuCountDisp= false)
end

overlayZTwindow(this, idx)

beautifyPlot()
% Adjust axes
xrange = opt.rangeWindow / (24*3600); % Convert seconds to days
xlim(xrange);

% Compute y-limit based on data inside the x-range only
inRange = t >= xrange(1) & t <= xrange(2);
ylim([0, max(freq(inRange))*1.1]);


title(sprintf('Rat %d - Zoomed IS count in time', opt.rat));
xlabel('Time (days)');
ylabel(sprintf('Normalized IS count per bin of %d s', opt.binCount));
legend();

end

%%        %        %         %        %    HELP FUNCTIONS    %         %        %        %         % 

function overlayZTwindow(this, idx, opt)
% overlayZTwindow - Affiche les phases nuit (light OFF) en zones grises
%
% Inputs :
%   this : EpyData object
%   idx  : index du rat (ligne dans refDatesTimestamps)
%   opt  : options avec ZTlightOn (default 8), ZTlightOff (default 20)

arguments
    this (1,1) EpyData
    idx (1,1) double 
    opt.ZTlightOn (1,1) double = 7
    opt.ZTlightOff (1,1) double = 19
end

tStart = this.refDatesTimestamps(idx,1);
tEnd   = this.refDatesTimestamps(idx,2);

% Liste des jours couverts
Days = dateshift(tStart, 'start', 'day') : days(1) : dateshift(tEnd, 'start', 'day');

hold on;
ylims = ylim(gca); % pour patch

for d = 1:length(Days)
    lightOff = Days(d) + hours(opt.ZTlightOff);      % ex: 20:00 le jour J
    lightOn  = Days(d) + days(1) + hours(opt.ZTlightOn); % ex: 08:00 le jour J+1

    % Tronquer aux bornes réelles de l'enregistrement
    t1 = max(lightOff, tStart);
    t2 = min(lightOn, tEnd);

    % Handle Label
    labelOpt = {'DisplayName','Light off'};
    if d > 1
        labelOpt = {'HandleVisibility','off'};
    end

    if t1 < t2
        % Conversion en secondes depuis tStart
        x1 = seconds(t1 - tStart)/(24*3600);
        x2 = seconds(t2 - tStart)/(24*3600);

        fill([x1 x2 x2 x1], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
            [0.2 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.2, labelOpt{:});
    end

    
end
end

