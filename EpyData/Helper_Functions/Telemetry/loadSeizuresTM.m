function seizures = loadSeizuresTM(session,verbose)
% Load seizure intervals from a telemetry session.
% 
% This function attempts to load a .seizures file containing seizure start/end times 
% (in milliseconds) for a given telemetry session. If the file is not found, a message 
% is displayed. If found, the data is parsed and converted into seconds.
%
% Syntax:
%   seizures = loadSeizuresTM(session)
%
% Input:
%   session - Full path to the telemetry session folder. Type: string or char.
%             Example: '/media/data-142/Rat986-20230502'
%   - verbose (opt)   - [boolean] (default = true). Display processing steps in command window.
%
% Output:
%   seizures - Nx2 array of seizure intervals [start_time, end_time] in seconds.
%
% Author: Corentin — 18/05/2025

arguments
    session (1,1) string
    verbose (1,1) {mustBeMember(verbose, [0,1])} = true
end

    f_eeg = 512;

    % Convert input to char and extract session ID
    session = char(session);
    [~, sessionID] = fileparts(session);
    cd(session);
    seizures = [];

    % Define expected filename
    seizure_file = [sessionID, '.seizures'];

    % Check if the .seizures file exists
    if isfile(seizure_file)
        % Read file contents and split into lines
        raw = fileread(seizure_file);
        lines = strsplit(strtrim(raw), newline);
        n = numel(lines);
        seizures = zeros(n, 2);

        % Parse each line as a pair of samples, then convert to seconds
        for i = 1:n
            values = sscanf(lines{i}, '%f,%f');
            if numel(values) == 2
                seizures(i, :) = values(:)' / f_eeg;
            else
                error("Malformed line in file: %s\n", lines{i});
            end
        end
        return
    else 
        if verbose
            fprintf('File %s not found.\n', seizure_file);
        end
    end
end