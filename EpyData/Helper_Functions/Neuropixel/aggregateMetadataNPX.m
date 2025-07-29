function [ratID_np, res_datetime, sessionPath] = aggregateMetadataNPX(path, ratID, date, realTime)
% getSessionMetadataNPX gathers metadata for all sessions of each Neuropixels rats. 
%
% USAGE:
%   (Expected called : with the direct ouput of getMetadataNPX method).
%   [ratID_np, res_datetime, sessionPath] = aggregateMetadataNPX(path, ratID, date, realTime)
%
% INPUTS:
%   - path: Cell array of path associated to the session.
%   - ratID: Cell array of rat IDs associated to the session. Each line is a session and must be linked to one rat. (e.g., {'596', '596', '597', ...}).
%   - date: Cell array of datetime objects associated to the session. Each line is a session and must be linked to one calendar date.
%   - realTime: Cell array of structures. Each element is a structure with a `.time` field (in seconds) giving informations about starting (ref to midnight) and duration of every segment of the session.
%               
%
% OUTPUT: vector of double, string or datetime with one line corresponding to one unique same rat (same order for every output).
%   - ratID_np : vector of double, ID of each different rat.
%   - res_datetime : Nx2 datetime vector. First column = Date and hour of the beggining of the first
%                    session of rats. Second column = NaT (future end time).
%   - sessionPath : cell of Nx1 string vector containing all the path sessions associated to one rat.

arguments
     path cell
     ratID cell 
     date cell 
     realTime cell
end

ratID = cell2mat(ratID); % Convert to char array if needed

% Identify unique rats
ratID_np = unique(ratID, 'stable');
n_rats = numel(ratID_np);
duration_days = zeros(n_rats, 1);
res_datetime = NaT(n_rats, 2);
sessionPath = {};

for i = 1:n_rats

    % Find all sessions for the current rat
    idx_sessions = find(ratID == ratID_np(i));

    sessionPath{i,1} = string(cellfun(@(x) x, path(idx_sessions)));

    % Retrieve corresponding dates and sort sessions chronologically
    dates_rat = date(idx_sessions);
    [~, sort_idx] = sort([dates_rat{:}]);
    idx_sessions = idx_sessions(sort_idx);

    % Get first session start datetime
    first_day = date{idx_sessions(1)};
    start_offset_s = realTime{idx_sessions(1)}.time(1); % s after midnight
    datetime_start = first_day + seconds(start_offset_s);

    res_datetime(i,1) = datetime_start;
end

end
%% Not exact endtime stored in the file RealTimeCat.evt (no gap taken in account) -> not used here
% % Get last session end datetime
% realTime = load(realtimecat.evt file)
% last_day = date{idx_sessions(end)};
% last_time_vector = realTime{idx_sessions(end)}.time;
% if last_time_vector(end-1) < 7*3600 % i. e. the beginning time of the last recording is before 7 am (arbitrary value) => the session has been restarted during the night of the next day
%     last_day = last_day + 1;
% end
% end_offset_s = last_time_vector(end-1) + last_time_vector(end) ; % s after midnight of the last day
% datetime_end = last_day + seconds(end_offset_s) ;
%
% res_datetime(i,2) = datetime_end;
% % Compute duration in days
% duration_days(i) = days(datetime_end - datetime_start);
