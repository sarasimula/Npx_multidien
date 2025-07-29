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
%   - RecDates (datetime): 1x2 vector - Date of the begging and ending recording (e.g., 2023-05-02).
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

%% Extraction of the precise begininng datetime

tm_recFolder = { % Based on the excel updated by Antoine : "/mnt/hubel-data-146/EdfData/bilan recordings OSI.xlsx"
    677,   ""; % No trace
    883,   "/mnt/hubel-data-146/EdfData/2022-12-Pilo/Alt157/";
    885,   "/mnt/hubel-data-146/EdfData/2022-12-Pilo/Alt157/";

    888,   "/mnt/hubel-data-146/EdfData/2022-12-Pilo/Alt160/";
    889,   "/mnt/hubel-data-146/EdfData/2022-12-Pilo/Alt160/";

    982,   "/mnt/hubel-data-146/EdfData/2023-05-Pilo/Alt157/";
    983,   "/mnt/hubel-data-146/EdfData/2023-05-Pilo/Alt157/";

    986,   "/mnt/hubel-data-146/EdfData/2023-05-Pilo/Alt160/";
    987,   "/mnt/hubel-data-146/EdfData/2023-05-Pilo/Alt160/";
    988,   "/mnt/hubel-data-146/EdfData/2023-05-Pilo/Alt160/";

    1065,  "/mnt/hubel-data-146/EdfData/2023-12-Pilo/Alt157/";
    1069,  "/mnt/hubel-data-146/EdfData/2023-12-Pilo/Alt157/";

    1070,  "/mnt/hubel-data-146/EdfData/2023-12-Pilo/Alt160/";    
    1074,  "/mnt/hubel-data-146/EdfData/2023-12-Pilo/Alt160/";

    1212,  "/mnt/hubel-data-146/EdfData/2024-06-Pilo/between_20240718_20240719/Alt157/";        
    1213,  "/mnt/hubel-data-146/EdfData/2024-06-Pilo/between_20240718_20240719/Alt157/";
    1214,  "/mnt/hubel-data-146/EdfData/2024-06-Pilo/between_20240718_20240719/Alt157/";

    1218,  "/mnt/hubel-data-146/EdfData/2024-06-Pilo/prior_20240718/Alt160/";
    1220,  "/mnt/hubel-data-146/EdfData/2024-06-Pilo/prior_20240718/Alt160/";
    1224,  "/mnt/hubel-data-146/EdfData/2024-06-Pilo/prior_20240718/Alt160/";

    1249,  "/mnt/hubel-data-146/EdfData/2024-10-Pilo/Partie1/Alt157/";
    1250,  "/mnt/hubel-data-146/EdfData/2024-10-Pilo/Partie1/Alt157/";

    1256,  "/mnt/hubel-data-146/EdfData/2024-10-Pilo/Partie1/Alt160/";    
    1257,  "/mnt/hubel-data-146/EdfData/2024-10-Pilo/Partie1/Alt160/";
};

% Find the good value in the cell and extract the good folder path
idx = find(sessionID==cell2mat(tm_recFolder(:,1))); 
folderPath = tm_recFolder{idx,2}; 

% Convert file name into datetime thanks to the unix code in the name of files and store in struct
if folderPath == "" % For rat 677 -> no trace of it
    startDatetime = recDates(1,1);
else % All the other rats
    files = dir(folderPath);
    fileNames = {files(~[files.isdir]).name};  % Exclusion of sub-folders

    % Alphanumerique sorting of the file names
    sortedNames = sort(fileNames);

    % Take the first name and extarct the unix to convert it in datetime
    startFile = char(sortedNames{1});
    unixTimestampSplitted = split(startFile(2:end), '.'); unixTimestamp = str2double(unixTimestampSplitted{1});
    startDatetime = datetime(unixTimestamp, 'convertfrom','posixtime','TimeZone','Europe/Madrid', 'Format','dd-MMM-yyyy HH:mm:ss.SSS');
    startDatetime.TimeZone = ''; 
end

% Sanity check between the datetime in the folder name and the one more
% precise
if  ymd(recDates(1)) ~= ymd(startDatetime)
    warning("CAREFUL FOR THE BEGINNING DATETIME")
end
recDates(1,1) = startDatetime;

%% Compute the cumulative timestamps for every segment and the end of the recording.
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
recDates(1,2) = recDates(1,1) + seconds(cumulativeTimes(end));

end
