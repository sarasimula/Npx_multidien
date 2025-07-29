function dataResc = rescaleNPX(session, data)
% rescaleNPX - For a folder, rescale Neuropixels session data on a continuous time axis in seconds
%
% USAGE:
%   dataResc = rescaleNPX(session, data)
%
% DESCRIPTION:
%   This function rescales timestamp data (e.g., ripple events or sleep intervals) recorded 
%   in multiple concatenated Neuropixels session segments (which constitue one folder) into a continuous time reference 
%   in seconds. It accounts for gaps between concatenated files, including crossing midnight,
%   and outputs corrected timestamps along with session start datetime and unrecorded intervals.
%
% INPUTS:
%   session : string or char
%       Full path to the session XML file (e.g., '/mnt/hubel-data-124/Rat596-20210618/Rat596-20210618.xml').
%
%   data : numeric matrix
%       Vector or matrix of timestamps in seconds (e.g., seizures, sleep intervals, IS timestamps ...).
%
% OUTPUT:
%   dataResc : struct with fields
%       .data     : matrix with rescaled timestamps (seconds), plus an extra column indicating day number (1 or 2).
%       .refStartAndStop    : 1x2 datetime vector of full datetime of the session start (from first segment)
%       i.e. reference for the beggining of timestamps and the full datetime of the ending time of the last segment.  
%       .unrecInt : Nx2 matrix of unrecorded intervals (in seconds) between concatenated segments.
%       will contain .ZT (if ZT activated) : double used as ZT. 
%
% NOTES:
%   - Requires '.RealTimeCat.evt' file in the same folder as session XML for event times and durations.
%   - Handles concatenated sessions spanning multiple days (e.g., >24h).
%   - Automatically adjusts small gaps (<60s) to 60s to correct possible timestamp shifts.



arguments
    session (1,1) string
    data double
end
dataresc = data;

%% Load Session
session = string(session);
[path, sessionID] = fileparts(session);

%% Load
% Load RealTimeCat event times (in seconds)
TimeConcat = LoadEvents(path + '/' + sessionID + '.RealTimeCat.evt'); % Load in Second !! (eventhought when opened displayed in ms)
N = length(TimeConcat.time);
Start = TimeConcat.time(1:2:N,1); % The first element of the time vector is a starting time of the first session of that day concatenated session, the second element is the duration of taht first session, the third is the starting time of the second session, the fourth is the suration of this second session and so on.
Duration = TimeConcat.time(2:2:N,1);

% Extract recording dates from descriptions
Dates = [];
for i=1:size(TimeConcat.description,1)/2
    Description = TimeConcat.description{2*i-1,1};
    dt = extractAfter(Description, '-data-');
    Dates = [Dates; datetime(dt(5:14),"InputFormat","uuuu-MM-dd")];
end

Z24 = 24*3600; % 24h in seconds

%% Compute gaps between concatenated segments (correct for crossing midnight)
% Gap between files : if 2 or more files
GapCon = [];

if size(Start,1) > 1
    for m = 1:length(Start)-1
        if day(Dates(m,1)) < day(Dates(m+1,1)) % If the recording spans between sessions on different days, we must account for crossing midnight (hence the +Z24 in the next line)
            GapCon(m) = Start(m+1) - (Start(m) + Duration(m)) + Z24; % Gap between files

        else % If the recording is within the same day, no additional correction is needed
            GapCon(m) = Start(m+1) - (Start(m) + Duration(m));
        end

        if GapCon(m)<60  % In case of a timestamp shifting (due to an unprecised sampling rate during recording), the gap between 2 sessions turns to be under 60 seconds (or even negative), we reset it to 60s (arbitrary value).
            GapCon(m) = 60;
        end
    end
else
    GapCon = nan; % No Gap if only one file
end

% Chech for physically impossible negative gap
if any(GapCon < 0) % With the if loop GapCon(m)<6à should not be usefull anymore
    disp(GapCon(GapCon < 0))
    error('Negative gaps detected (may due to incoherent starting point) for session %s:', sessionID);
end

%% Correction
% 1 - Correction by the concatenation Gap (GapCon)

if ~isnan(GapCon) % Only if concatenated recordings because sometimes only one recording
    Ttransition = cell(1, size(dataresc,2));  % Initialisation
        for col = 1:size(dataresc,2)
        transitions = [];
            % find T transition between recording
            for m = 1:size(GapCon,2) % if more than 2 recordings
                idx = find(dataresc(:,col) > sum(Duration(1:m)), 1, 'first');
                if ~isempty(idx)
                    transitions(end+1) = idx;  % Add only if found
                end
            end
        Ttransition{col} = transitions;
        end
        
        if size(Ttransition,2)>1 % In case of a data with 2 vectors (intervals of sleep for example)
            duration = 0;
            col1 = dataresc(:,1);
            col2 = dataresc(:,2);
            for i = 1:numel(Ttransition{1})
                duration = duration + Duration(i);
                if  Ttransition{1}(i) > Ttransition{2}(i) % if an interval is overlapping 2 segments (because the event detections have been done on the concatenated segments so it's possible and it's happenning for sleep intervals).
                    % insertion of the end of the first segment in both columns at their right idx
                    col1 = sortrows([col1;duration]);
                    col2 = sortrows([col2;duration]); 
                    Ttransition{2}(i) = Ttransition{2}(i)+1;
                    
                    if i < numel(Ttransition{1}) % Because we added a new row -> readjust the previous row index values.
                        Ttransition{1}(i+1:end) = Ttransition{1}(i+1:end) +1;
                        Ttransition{2}(i+1:end) = Ttransition{2}(i+1:end) +1;
                    end
                end
            end
            dataresc = [col1,col2];
        end
        


    % Correction
    for col = 1:size(dataresc,2)
        transitions_l = Ttransition{col};
        for k = 1:length(transitions_l)
            idx = transitions_l(k);
            if ~isnan(idx)
                dataresc(idx:end, col) = dataresc(idx:end, col) + GapCon(k);
            end
        end
    end

end


% Compute the ending time of the session 
endTime = Dates(1) + seconds(Start(1)) + seconds(Duration(1));
if ~isnan(GapCon)
    for segment = 1:size(GapCon,2)
        endTime = endTime + seconds(GapCon(segment))+seconds(Duration(segment+1));
    end
end

% Save the unrecorded time intervals in the session
unrecorded_intervals = [];
duration = 0;
for col = 1:length(Start)-1
    duration = duration + Duration(col);
    unrecorded_intervals(col,1) = duration;
    duration = duration + GapCon(col);
    unrecorded_intervals(col,2) = duration;
end
% [path, sessionID, ~] = fileparts(session);
% output_filename = fullfile(path, [sessionID '_unrecorded_intervals.txt']);
% writematrix(unrecorded_intervals, output_filename, 'Delimiter', '\t');




%% Output
dataResc.data = dataresc;
dataResc.refStartAndStop = [Dates(1) + seconds(Start(1)), endTime] ;
dataResc.unrecInt = unrecorded_intervals;

end

