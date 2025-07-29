function [path, ratID, date, realTime] = getMetadataNPX(session)
%GETSESSIONNPX Returns informations f one session for NPX data.
%
% USAGE:
%   [path, ratID, date, realTime] = getMetadataNPX(session)
%
% INPUT:
%   - session (char/string): Path to a session file (Example of session '/mnt/hubel-data-124/Rat596-20210618/Rat596-20210618.xml' ).
%
% OUTPUT:
%   - path (string): Path to the session directory.
%   - ratID (double) : ID of the rat
%   - date (datetime) : calendar date time of the session
%   - realTime (struct) : tracke the beggining with fields
%           - .time (vector of double) : The first element (and globaly the even index element) of the time vector is a starting time of the first (or second, third ...) segment of that day concatenated session
%                     The second element (and globaly the odd index element) is the duration of that first (or second, third ...) segment.
%           - .description : validation of the vector time (unused)


arguments
     session (1,1) string
end


    %% Extract Rat ID and Session Date
    session = string(session);
    [path, sessionID] = fileparts(session);

    % Example: sessionID = 'Rat0885-20221216'
    ratID = str2double(cell2mat(extractBetween(sessionID, 'Rat', '-')));
    dt = extractAfter(sessionID, '-');
    date = datetime(dt,"InputFormat","uuuuMMdd");

        %% Load Real Time Events
    % Load synchronization timestamps from the .RealTimeCat.evt file
    try 
        realTime = LoadEvents(strcat(path,'/',sessionID,'.RealTimeCat.evt'));
    catch
        error("No .RealTimeCat.evt for %s", session)
    end

end