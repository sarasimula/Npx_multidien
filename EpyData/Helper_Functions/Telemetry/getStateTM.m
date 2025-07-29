function [sleepInt, swsInt, wakeInt] = getStateTM(session,tsMicroAwake, verbose)
% Extracts and concatenates sleep intervals (from .sws and .rem file) and deduces wake intervals from a telemetry rat.
%
% Syntax:
%   [sleepInt, wakeInt] = getStateTM(session,micro_awakeness_treshold)
%
% Input:
%   session - Full path to a telemetry session folder. Type: string or char.
%             Example: '/media/data-142/Rat986-20230502'
%   tsMicroAwake (opt) - double. Treshold in second, under the value it erases micro awakeness intervals and convert it into sleep time.
%   verbose (opt)   - [boolean] (default = true). Display processing steps in command window.
%
% Output:
%   sleepInt - Nx2 array of sleep intervals (in seconds) computed over the session.
%   swsInt   - Nx2 array of sws intervals (in seconds) computed over the session.
%   wakeInt  - Nx2 array of wake intervals (in seconds) deduced over the session.
%
% Author: Corentin — 21/05/2025 (last edit : ) 

arguments
    session (1,1) string
    tsMicroAwake (1,1) double = 0; % Treshold in second, under the value it erases micro awakeness intervals and convert it into sleep time.
    verbose (1,1) {mustBeMember(verbose, [0,1])} = true
end

% Get session ID and JSON metadata
session = string(session);
[~, sessionID] = fileparts(session);
cd(session)

% Test if presence of sleep
[sleepInt, swsInt, wakeInt] = deal([]);
if isempty(dir('*.sws')) && isempty(dir('*.rem'))
    if verbose
        fprintf("No sleep file found for this rat : %s\n", char(session));
    end
    return;
end

% Load JSON metadata
jsonPath = session + '/' + sessionID + '_sessions_total_points.json';
fid = fopen(jsonPath);
raw = fread(fid,  inf, 'uint8=>char')';
fclose(fid);
fn = jsondecode(char(raw));

% Extract duration of each segment
points = cell2mat(struct2cell(fn));

% Compute segment start times in seconds (ccumlative times)
freq_eeg = 512; % Hz
cumulativeTimes = cumsum([0; points(1:end)]) / freq_eeg; % Each elements i give the number of seconds spent from the beggining to the session i. The last element of cumulative_times give the end of recording.

% Extract sleep periods
for i = 1:numel(points)
    if isempty(dir("*"+ num2str(i-1) +".sws")) && isempty(dir("*"+ num2str(i-1) +".rem")) 
        if verbose
            disp("No sleep data for segment " + num2str(i-1) + " of session " + sessionID);
        end
    else
        sws = [];
        rem = [];
        try 
            sws = dlmread(fullfile(session, sessionID + "_" + num2str(i-1) + ".sws"));
        catch
            if verbose
                disp("No .sws (but presence of .rem) file for segment " + num2str(i-1) + " of session " + sessionID);
            end
        end

        try 
            rem = dlmread(fullfile(session, sessionID + "_" + num2str(i-1) + ".rem"));
        catch
            if verbose
                disp("Warning : no .rem (but presence of .sws) file for segment " + num2str(i-1) + " of session " + sessionID);
            end
        end
        
        if ~isempty(rem) && ~isempty(sws)
            sleep = sortrows([sws; rem]);
        elseif isempty(rem)
            sleep = sws;
        elseif isempty(sws)
            sleep = rem;
        end

        % Compute the wake -> for telemtry ther is only partition in two
        % states : sleep and wake ther is no unrecorded intervals.
        if ~isempty(sleep) && ~isempty(sws) % We need to have at least one sleep intervals during a segment. It should always be the case if we are at this stage (because we loaded at least .rem or .sws), but it appears some segments have a .sws or a .rem but the file is empty. Example : no .rem (but presence of .sws) file for segment 8 of session Rat1212-20240701.
            % To compute wake, we want at least a functionnal sws file (more intervals), just a rem intervals is not a good approximation of the amount of sleep.

            % Supress micro awakeness
            sleep = cleaned_intervals(sleep,tsMicroAwake); % Erase micro awakenes (< treshold) and convert it into sleep time.

            windowTimestamps = [cumulativeTimes(i) , cumulativeTimes(i+1)];
            wake = zeros(size(sleep,1)+1,2);
            wake(1,1) = windowTimestamps(1);
            wake(end,2) = windowTimestamps(2);
            wake(2:end,1) = sleep(:,2) + cumulativeTimes(i);
            wake(1:end-1,2) = sleep(:,1) + cumulativeTimes(i) ;
            wakeInt = [wakeInt; wake];

            % Debugging
            % checkIntervals(sessionID + sprintf('_segment%d',i-1), sleep, 1*60*60+30*60,20)
            % checkIntervals(sessionID + sprintf('_segment%d',i-1), wake, 4*60*60+30*60,20)

        end

        % Store the state for this segment
        sleepInt = [sleepInt; sleep + cumulativeTimes(i)];
        swsInt = [swsInt; sws + cumulativeTimes(i)];
       

    end
end



end
