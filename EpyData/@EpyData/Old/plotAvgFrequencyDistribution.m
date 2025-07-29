function xtick = plotAvgFrequencyDistribution(this, opt)
%PLOTAVG DISTRIBUTION Plot a boxplot of average IS or Seizure rate per rat.
%
% Syntax:
%   obj.plotAvgDistribution(opt)
%
% Arguments:
%   opt.nature (1,1) string {mustBeMember(opt.nature, ["IS", "Seizures"])} = "IS"
%   opt.state  (1,1) string {mustBeMember(opt.state, ["all", "sleep", "wake"])} = "all"
%   opt.figRatio (1,1) double = 1.5

arguments
    this (1,1) EpyData
    opt.rat(1,:) double = this.ratID  % optionnal, by défault all the rats
    opt.nature (1,1) string {mustBeMember(opt.nature, ["IS", "seizures"])} = "IS"
    opt.state (1,1) string {mustBeMember(opt.state, ["all", "sleep", "wake"])} = "all" % By default, all IS
    opt.binCount (1,1) double = 3600  % Bin size in seconds for the count
    opt.range (1,1) double {mustBeNonnegative} = 0
    opt.color (3,1) double = [0.3, 0.3, 0]
    opt.xtick (1,1) double = 1  % Position of the violon on x axis
    opt.figRatio (1,1) double = 1.5
end

%% === EXTRACTION DES DONNÉES ===
% Extract the avg_vector from `this` and from options `opt`
count = this.getFrequency("nature",opt.nature, "rat", opt.rat,"state",opt.state,"binCount",opt.binCount, "range",opt.range);
avg_vector = count.avg;  % Ex: avg IS/sec per rat, or seizures/sec, depending on opt

%% === PLOT UNIQUE : BOXPLOT ===
% Dimensions
height_cm = 10;
width_cm = opt.figRatio * height_cm;

% Plot
xtick = opt.xtick;
daviolinplot(avg_vector, 'x', opt.xtick, ...
    'colors', opt.color, 'boxcolors', 'k', ...
    'outliers', 0, 'box', 0, 'boxwidth', 0.8, ...
    'scatter', 2, 'scattersize', 14, 'jitter', 1, ...
    'scattercolors', opt.color, 'bins', 12);


ylabel(sprintf('Average %s rate (Hz)', opt.nature), 'Interpreter', 'none');
grid on;

set(gca, 'FontSize', 10);


end
