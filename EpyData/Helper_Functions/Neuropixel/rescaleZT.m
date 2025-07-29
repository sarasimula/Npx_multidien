function  dataZT = rescaleZT(data, ZT)
% rescaleZT - Rescale timestamps relative to Zeitgeber Time (ZT) reference
%
% USAGE:
%   dataZT = rescaleZT(data, ZT)
%
% DESCRIPTION:
%   Rescales timestamps (in seconds) relative to a Zeitgeber Time (ZT) zero reference,
%   which is the beginning of the light phase, expressed in hours (e.g., ZT = 7 means light phase starts at 7:00).
%   Adjusts timestamps so that times before ZT belong to the previous day, 
%   thus aligning data on the ZT clock.
%
% INPUTS:
%   data : struct
%       .data     : matrix of timestamps in seconds (output from rescaleNPX)
%       .refStartAndStop    : 1x2 datetime vector [full datetime of the data start, full datetime of the data ending]
%       .unrecInt : unrecorded intervals (output from RescaleNPX)
%
%   ZT : numeric scalar
%       Zeitgeber Time zero in hours (e.g., 7 means 7 AM)
%
% OUTPUT:
%   dataZT : struct
%       .data     : matrix of timestamps rescaled on Zeitgeber time (seconds)
%       .refStartAndStop    : 1x2 datetime vector [new full datetime of the data start based on the ZT, same full datetime of the data ending]
%       .unrecInt : rescaled unrecorded intervals
%       .ZT       : double used as ZT.
%
% NOTES:
%   - This function expects 'data' to be the output struct from rescaleNPX or
%     equivalent struct containing 'data' and 'refStartAndStop'.
%   - Adjusts the reference start date/time to align with ZT zero.

arguments
    data struct
    ZT (1,1) double
end

% Get date of session start, reset to midnight
dataZT.ZT = ZT;
startSession = data.refStartAndStop(1);
date = dateshift(startSession,'start','day');

%If session start hour is before ZT, shift date backward by one day
if hours(startSession - date) < ZT
    date = date - days(1); % previous day for alignment
end

% Define new reference start datetime at ZT hour
dataZT.refStartAndStop(1) = date + hours(ZT);
dataZT.refStartAndStop(2) = data.refStartAndStop(2); % No change of the ending

% Calculate the offset in seconds between original start and new ZT reference
secGapFromZT = seconds(startSession - dataZT.refStartAndStop(1));

% Rescale timestamps accordingly
dataZT.data = data.data + secGapFromZT;
dataZT.unrecInt = data.unrecInt + secGapFromZT;

end

%% Identify if timestamps cross into a second day (>24h) and add day number column
% 
% Tday = [];
% 
% data = dataZT.data;
% 
% if isdmatrix(data)
%     for n = 1:size(data,2)
%         y = find(data(:,n) > Z24,1,'first'); % find T if longer than 24h (following day)
%         Tday = [Tday,y];
%     end
% else
%     Tday = find(data > Z24,1,'first');
% end
% 
% data = [data ones(size(data,1),1)];
% 
% if ~isempty(Tday)
%     data(Tday(1):end,size(data,2)) = 2;
% end
% 
% dataZT.data = data;