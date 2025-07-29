function h_all = getHistogramIS(this, range, opt) 
% This function computes a frquency of events based on a 
% defined time bin step and a max limit (range). The output represents the count of events 
% per time bin, with zero values replaced by NaN for improved visualization.
%
% Syntax :  h_all = get_histogram_PSD_Corentin(all_IS, histogram_step,range) 
%
% Inputs :
%   all_IS - Timestamps of events (in seconds). Type: Vector.
%   histogram_step  - Time bin size for the histogram (in seconds). Type: Scalar.
%   range  - Upper time limit (in seconds) defining the maximum histogram bin value. Type: Scalar. 
%
% Outputs :
%   h_all  - Histogram of event counts per time bin. Zeros are replaced by NaN. Type: Vector.
%
% Authors : Corentin (original: 13/03/2025 . Last edit: 13/03/2025 .),  
% Contributors : 

% DEBUG/
% all_IS = randi([0, 86400], 1, 100); % Example: 100 random IS timestamps within a day (in seconds)
% step = 3600; % Example: 1-hour bin size
% h_all = get_histogram_PSD_Corentin(all_IS, step);

arguments
  this (1,1) EpyData
  range (1,1) double {mustBePositive} = max(this.ISallTime)
  opt.step (1,1) double {mustBePositive} = 3600
  opt.state (1,1) string {mustBeMember(opt.state, ["all", "sleep", "wake"])} = "all"
end

% Histogram Analysis of Interictal Spikes (IS)
edges = 0:opt.step:floor(range); % Define histogram edges, range defines the top limit for time bins

% Compute IS counts (per step)
if opt.state == "all"
    h_all = histcounts(this.ISallTime, edges); % Compute histogram
elseif opt.state == "wake"
    h_all = histcounts(this.ISwakeTime, edges);
elseif opt.state == "sleep"
    h_all = histcounts(this.ISsleepTime, edges);
end

h_all(h_all == 0) = NaN; % Replace zeros with NaN for better visualization
end