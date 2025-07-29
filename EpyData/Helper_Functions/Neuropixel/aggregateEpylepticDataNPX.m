function [ISallTime, seizuresTime, sleepInt, swsInt, unrecInt, wakeInt, refDatesSessionTimestamps, refDatesTimestamp, ZT, UnitsallTime] = aggregateEpylepticDataNPX(ratID, IS_rescale, seizures_rescale, sleep_rescale, sws_rescale, wake_rescale, units_rescale)
% From raw structures, aggregate informations (IS, seizures, sleep, unrec, wake)
%
% INPUTS:
%   - ratID: Cell array containing the rat ID for each session.
%   - date: Cell array with MATLAB datetime objects for each session.
%   - IS_rescale, seizures_rescale, sleep_rescale, wake_rescale (output of loadEpilepticEventsNPX) : All cell arrays of Nx1 structures with the same fields :
%       .data     : matrix with rescaled timestamps (seconds)
%       .unrecInt : Nx2 matrix of unrecorded intervals (in seconds) between concatenated segments.
%       .refStartAndStop    : 1x2  full datetime vector[start, stop] i.e. references for the beggining and ending of all timestamps of the struct.
%       may contain .ZT : double used as ZT.
%
% OUTPUTS:
%   - ISallTime: Cell array of all IS timestamps (in seconds) per rat
%   - seizuresTime: Cell array of Nx2 vector of seizures timestamps (in seconds) intervals per rat
%   - sleepInt: Cell array of Nx2 vector of sleep timestamps (in seconds) intervals per rat
%   - swsInt : Cell array of Nx2 vector of slow waves sleep timestamps (in seconds) intervals per rat
%   - unrecInt: Cell array of Nx2 vector of unrec timestamps (in seconds) intervals per rat
%   - wakeInt: Cell array of Nx2 vector of wake timestamps (in seconds) intervals per rat
%   - refDatesSessionTimestamp: Cell array of stuctures with two fileds :
%           .date : Nx2 full datetime vector for beggining and ending of each session
%           .cumulative : Nx2 double vector with cumulative timestamp in seconds from reference time given by refDatesTimestamp(:,1)
%   - refDatesTimestamp: Nx2 full datetime vector with, per lines, dates of the begging and ending recording (format, 2023-05-02 10:52:03) of the entire recording.
%   - ZT: [int] Zeitgeber Time (default ZT=NaN : no rescaling to ZeitGeber Time) with which all the timestamps have been rescaled on.

arguments
    ratID cell 
    IS_rescale cell
    seizures_rescale cell
    sleep_rescale cell
    sws_rescale cell
    wake_rescale cell
    units_rescale cell = {}
end  

ratID = cell2mat(ratID);
ratID_np = unique(ratID, 'stable'); % Unique rat IDs, preserving order
n_rats = numel(ratID_np); % Number of distinct rats

% Initialize output cell arrays
[ISallTime, seizuresTime, sleepInt, swsInt, wakeInt, unrecInt, refDatesSessionTimestamps, UnitsallTime] = deal(cell(n_rats,1));
refDatesTimestamp = NaT(n_rats,2);

for i = 1:n_rats

    % Initialize arrays for one rat across all its recording days
    [ISallTime_rat, seizuresTime_rat, sleepInt_rat, swsInt_rat, wakeInt_rat, unrecInt_rat, cumuStartFromFirstSession, cumuEndFromFirstSession, UnitsallTime_rat] = deal([]);
   
    Idx_days_rat_i = find(ratID == ratID_np(i)); % This vector’s length gives the number of days for rat i. Its values indicate the indices of the sessions for that rat.
    refDatesSessionTimestamps_rat = struct('date', NaT(numel(Idx_days_rat_i),2), 'cumulative', NaN(numel(Idx_days_rat_i),2));

    % initialize count of clusters ID to differentiate them day by day
    unitscount=0;

    % [refDatesSessionTimestamps_rat.date, refDatesSessionTimestamps_rat.cumulative] = deal(NaT(numel(Idx_days_rat_i),2));

    % Update the cumulative time for each day composing the sessions of rat i
    for n = 1:length(Idx_days_rat_i)
        
        % Sanity check of the ref time for all the 
        % TBD

        % The only type of file we are sure it existed for all session is IS (seizures and sleep are not always computed)
        cumuStartFromFirstSession(n) = seconds(IS_rescale{Idx_days_rat_i(n)}.refStartAndStop(1) - IS_rescale{Idx_days_rat_i(1)}.refStartAndStop(1));

        %% IS 

        % All IS
        is_session = IS_rescale{Idx_days_rat_i(n)}.data;
        ISallTime_rat = [ISallTime_rat; is_session + cumuStartFromFirstSession(n)];

        %% Units
        if ~isempty(units_rescale)
            units_session =  units_rescale{Idx_days_rat_i(n)}.data(:,1);
            cluname = units_rescale{Idx_days_rat_i(n)}.data(:,2);
            cluname = cluname + unitscount;
            unitscount = unitscount + max(cluname);
            UnitsallTime_rat = [UnitsallTime_rat; [units_session + cumuStartFromFirstSession(n), cluname]];
            clearvars units_session cluname
        else, UnitsallTime_rat = NaN; 
        end
        

        %% Sleep & wake
        % Sanity check
        if ~ isequal(fieldnames(sleep_rescale{Idx_days_rat_i(n)}), fieldnames(wake_rescale{Idx_days_rat_i(n)}))
            error('Wake intervals NOT WELL IMPLEMENTED -> ISSUE')
        end

        % Sleep & wake intervals if existing
        if ~isempty(fieldnames(sleep_rescale{Idx_days_rat_i(n)})) ... % If there are intervals of sleep for the session (so intervals of wake)
                && ~isempty(fieldnames(wake_rescale{Idx_days_rat_i(n)})) % Not useful for computation reason but better for understanding

            sleep_session = sleep_rescale{Idx_days_rat_i(n)}.data; % Matrix containing start and end timestamps of each sleep episode
            sleepInt_rat = [sleepInt_rat; sleep_session + cumuStartFromFirstSession(n)];

            wake_session = wake_rescale{Idx_days_rat_i(n)}.data;
            wakeInt_rat = [wakeInt_rat; wake_session + cumuStartFromFirstSession(n)];

        end

        % SWS intervals
        if ~isempty(fieldnames(sws_rescale{Idx_days_rat_i(n)})) % If there are intervals of sleep for the session (so intervals of wake)
            sws_session = sws_rescale{Idx_days_rat_i(n)}.data; % Matrix containing start and end timestamps of each sleep episode
            swsInt_rat = [swsInt_rat; sws_session + cumuStartFromFirstSession(n)];
        end

        %% Seizures
        % Seizure extraction for day n
        if ~isempty(fieldnames(seizures_rescale{Idx_days_rat_i(n)}))
            seizures_session = seizures_rescale{Idx_days_rat_i(n)}.data;
            seizuresTime_rat = [seizuresTime_rat; seizures_session + cumuStartFromFirstSession(n)];
        end

        %% UnrecInt
        % Unrecorded intervals beetween different sessions
        cumuEndFromFirstSession(n) =  seconds(IS_rescale{Idx_days_rat_i(n)}.refStartAndStop(2) - IS_rescale{Idx_days_rat_i(1)}.refStartAndStop(1) );
        if n>1 % Exclude the first session (no gap between an other interval before cause it's the first)
            unrecInt_rat = [unrecInt_rat; cumuEndFromFirstSession(n-1) cumuStartFromFirstSession(n)];
        end

        % Unrecorded intervals from the session n due do to different segments
        unrecInt_session = IS_rescale{(Idx_days_rat_i(n))}.unrecInt;
        unrecInt_rat = [unrecInt_rat; unrecInt_session + cumuStartFromFirstSession(n)];

        %% Save the beggining and ending timestamps (real date and cumulative time) of every session
        refDatesSessionTimestamps_rat.date(n,:)  = [IS_rescale{Idx_days_rat_i(n)}.refStartAndStop(1) IS_rescale{Idx_days_rat_i(n)}.refStartAndStop(2)];
        refDatesSessionTimestamps_rat.cumulative(n,:)  = [cumuStartFromFirstSession(n) cumuEndFromFirstSession(n)];

    end

    
    %% Output storage
    ISallTime{i,1} = ISallTime_rat;
    seizuresTime{i,1} = seizuresTime_rat;
    sleepInt{i,1} = sleepInt_rat;
    swsInt{i,1} = swsInt_rat;
    wakeInt{i,1} = wakeInt_rat;
    unrecInt{i,1} = unrecInt_rat;
    refDatesSessionTimestamps{i,1} = refDatesSessionTimestamps_rat;
    refDatesTimestamp(i,:) = [IS_rescale{Idx_days_rat_i(1)}.refStartAndStop(1), IS_rescale{Idx_days_rat_i(numel(Idx_days_rat_i))}.refStartAndStop(2)]; % SUPPOSE A CROISSANT SORT OF THE BATCH : Timestamps of the reference Dates for the rat n
    ZT = IS_rescale{Idx_days_rat_i(1)}.ZT;
    UnitsallTime{i,1} = UnitsallTime_rat;
end

end

    % %% Compute wake intervals
    % % Get the full session window : [start of the first segment of the first day, end of the last segment of the ultimate day]
    % refDatesTimestamp_rat = [IS_rescale{Idx_days_rat_i(1)}.refStartAndStop(1), IS_rescale{Idx_days_rat_i(numel(Idx_days_rat_i))}.refStartAndStop(2)];
    % windowTimestamps = [0 seconds(refDatesTimestamp_rat(2)-refDatesTimestamp_rat(1))];
    % stateVect = sortrows([sleepInt_rat;unrecInt_rat]);
    % wakeInt_rat = zeros(size(stateVect,1)+1,2);
    % wakeInt_rat(1,1) = windowTimestamps(1);
    % wakeInt_rat(end,2) = windowTimestamps(2);
    % wakeInt_rat(2:end,1) = stateVect(:,2);
    % wakeInt_rat(1:end-1,2) = stateVect(:,1);
