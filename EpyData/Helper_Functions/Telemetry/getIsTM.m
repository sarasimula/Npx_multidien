function [all_IS] = getIsTM(session, verbose)
%GETISTM Extracts interictal spike (IS) timestamps from a telemetry session.
%
% This function processes IS event files from telemetry rat recordings and returns
% all IS timestamps.
%
% Syntax:
%   [all_IS, all_sleep_IS, all_wake_IS] = getIsTM(session, wakeAndSleepIsSeperated)
%
% Inputs:
%   - session (string): Full path to the session folder (e.g., '/media/data/Rat677-20230502').
%   - verbose (opt)   - [boolean] (default = true). Display processing steps in command window.
%
% Outputs:
%   - all_IS        (vector): IS timestamps (in seconds) across the session.
%
% Authors : Corentin (original: 13/03/2025 . Last edit: 02/04/2025 .),  
% Contributors : 

% DEBUG/
% session = '/mnt/hubel-data-146/Rat0883-20221216';

arguments
    session (1,1) string
    verbose (1,1) {mustBeMember(verbose, [0,1])} = true
end

% Extract the session ID and open the session path
session = string(session);
[~, sessionID] = fileparts(session);
cd(session)

% Read the JSON file containing the number of time unit for each session (from 0 to X sessions)
fname = session + "/" + sessionID + "_sessions_total_points.json"; 
fid = fopen(fname); 
raw = fread(fid,  inf, 'uint8=>char')';
fclose(fid);
fn = jsondecode(char(raw));

% Extract duration of each segment
points = cell2mat(struct2cell(fn));

% Compute cumulative timeline for IS events
freq_eeg = 512; % EEG sampling frequency (Hz, samples per second)
cumulative_times = cumsum([0; points(1:end-1)]) / freq_eeg; % Each elements i give the number of seconds spent from the beggining to the session i. The last element of cumulative_times give the end of recording.

%% Load and Process Interictal Spike (IS) Data
% Initialize storage arrays for IS occurrences
all_IS = [];       % Timestamps (in seconds) of all IS events 
 
for i = 1:numel(points)-1
    try
        IS = dlmread(sessionID + "_" + num2str(i-1) + ".IS_Q60");
        all_IS = [all_IS; IS + cumulative_times(i)]; % Timestamps (in seconds) of all IS events
    catch
        if verbose
            disp("IS file is missing for segment " + num2str(i-1) + " for session " + sessionID);
        end
    end
end
end


