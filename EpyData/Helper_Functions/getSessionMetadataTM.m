function [sessionPath, sessionID, recDates, cumulativeTimes] = getSessionMetadataTM(session)
%GETSESSIONMETADATA Extracts metadata and recording duration for a telemetry rat session.
%
% This function is designed for telemetry recordings of rats. It extracts metadata
% such as the session path, numeric session ID, and recording date, and computes
% the total duration of the recording in days based on a JSON file containing
% the number of data points recorded per day.
%
% INPUT:
%   - session (char/string): Full path to the session folder, e.g.,
%       '/mnt/hubel-data-146/Rat677-20230502'
%
% OUTPUT:
%   - sessionPath (string): Full session path.
%   - sessionID (double): Numeric ID of the rat (e.g., 677).
%   - RecDates (datetime): 1x2 vector - Date of the begging and ending recording (e.g., 2023-05-02 10:52:03).
%   - CumulativeTimes (double) : Nx1 - Cumulative number of seconds spent
%   everyday with the reference to the beggining of the session. The last
%   element is the ending timestamp of the recording.
%
% NOTE:
%   This function requires the presence of the JSON file
%   '<sessionID>_sessions_total_points.json' in the session folder,
%   which must contain the number of recorded points per day.

arguments
    session (1,1) string
end

% Constants
freq_eeg = 512; % Hz
sec_per_day = 86400; % seconds in one day

recDates = NaT(1,2);

% Normalize session path
sessionPath = string(session);
[~, sessionData] = fileparts(sessionPath);

% Extract numeric session ID and recording date using regex
tokens = regexp(sessionData, "Rat(\d+)-(\d{8})", "tokens");
tokens = tokens{1}; % extract match
sessionID = str2double(tokens{1});
recDates(1,1)   = datetime(tokens{2}, 'InputFormat', 'yyyyMMdd');

% Construct full path to the JSON file
json_filename = sessionPath + '/' + sessionData + '_sessions_total_points.json';

% Read JSON content
if ~isfile(json_filename)
    error('JSON file not found: %s', json_filename);
end
fid = fopen(json_filename);
raw = fread(fid, inf, 'uint8=>char')';
fclose(fid);
fn = jsondecode(char(raw));
points = cell2mat(struct2cell(fn));

% Compute cumlative recorded time in seconds based on the sampling rate and the sampling points per day
cumulativeTimes = cumsum([0; points(1:end)]) / freq_eeg; % Each elements i give the number of seconds spent from the beggining to the session i. The last element of cumulative_times give the end of recording.

% Convert to days
nRecDays = cumulativeTimes(end) / sec_per_day;
recDates(1,2) = recDates(1,1) + days(nRecDays);

end
