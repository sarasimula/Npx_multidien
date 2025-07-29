function [IS, seizures, sleep, sws, wake, ISType, ISWaveforms, singleunits] = getEpilepticEventsNPX(session, tsMicroAwake, ZT, verbose, units)
% getEpilepticEventsNPX  Load epileptic event data from Neuropixels session.
%
% USAGE:
%   [IS, seizures, sleep, ISType, ISWaveforms] = getEpilepticEventsNPX(session, ZT, verbose)
%
% This function extracts various time-aligned event markers (IS, seizures, sws) and
% real-time reference from a session folder, and rescales them to Zeitgeber Time.
%
% Syntax : [ratID, date,realTime,IS_rescale,seizures,sws,waveforms,is_type] = LoadEpilepticEvents_Corentin(session)
%
% INPUT:
%   session             - [char or string] Full path to the session folder, or session file path.
%   ZT (opt)            - [int] Zeitgeber Time (default ZT=NaN : no rescaling to ZeitGeber Time) to rescale data on.
%   verbose (opt)       - [boolean] (default = true). Display processing steps in command window.
%   tsMicroAwake (opt)  - [double]. Treshold in second, under the value it erases micro awakeness intervals and convert it into sleep time.
%
% OUTPUT:
%   IS - [struct] Aligned timestamps (in seconds) with file time gaps (.realTimeCat.evt) of Interictal Spikes (IS), from .cortex_IS_corr.
%                        The time reference for the rescaling is by default the beginninig of the recording or the ZeitGeber time if given.
%                        
%   seizures   - [struct] Aligned timestamps (in seconds) with file time gaps (.realTimeCat.evt) of epileptic seizures intervals from .seizures.
%                        The time reference for the rescaling is by default the beginninig of the recording or the ZeitGeber time if given.
%
%   sleep      - [struct] Aligned timestamps (in seconds) with file time gaps (.realTimeCat.evt) of sleep intervals (rem + sws) from .rem & .sws.
%                        The time reference for the rescaling is by default the beginninig of the recording or the ZeitGeber time if given.
%
%   sws        - [struct] Aligned timestamps (in seconds) with file time gaps (.realTimeCat.evt) of sleep intervals from .sws.
%                        The time reference for the rescaling is by default the beginninig of the recording or the ZeitGeber time if given.
%
%   wake       - [struct] Aligned timestamps (in seconds) with file time gaps (.realTimeCat.evt) of wake intervals (thanks unrecorded intervals and sleep intervals).
%                        The time reference for the rescaling is by default the beginninig of the recording or the ZeitGeber time if given.
%
%   NB : All the structures have the same fields :
%       .data     : matrix with rescaled timestamps (seconds)
%       .unrecInt : Nx2 matrix of unrecorded intervals (in seconds) between concatenated segments.
%       .refStartAndStop    : 1x2 datetime vector of full datetime of the session start (from first segment)
%       i.e. reference for the beggining of timestamps and the full datetime of the ending time of the last segment.  
%
%
%   TO BE ADJUST (Sara) : 
%   ISType    - [cell] IS type (0=positive IS, 1= neg IS with afterwave,2=neg IS without afterwave). Time in seconds.
%   waveforms  - [cell] IS waveforms (over 1250 samples)
%
% Edited by Corentin & Sara
% on 10/04/25: added waveforms and is_type output; 
% on 15/04: added rem file to sleep intervals (before: only sws)


arguments
    session (1,1) string
    tsMicroAwake (1,1) {mustBeNonnegative(tsMicroAwake)} = 0
    ZT double = NaN % ZeitGebber for NPX
    verbose (1,1) {mustBeMember(verbose, [0,1])} = true
    units (1,1)  {mustBeMember(units, [0,1])} = false
end

    % Get path and session ID
    session = string(session); 
    [path, sessionID] = fileparts(session);

    %% Load units
    if units
        SetCurrentSession(session)
         
        clunums = GetUnits; % 1st column= channel num, 2nd= cluster num (=unit)
        s = GetSpikeTimes('output','numbered');
    
        singleunits = rescaleNPX(session, s(:,1));
        singleunits.data(:,2) = s(:,2);

    else 
        singleunits = NaN;
    end
   
    %% Load Interictal Spikes (IS) and rescale it inside the session
    % Load raw IS timestamps (in seconds) from .cortex_IS_corr
    IS = dlmread(strcat(path,'/',sessionID,'.cortex_IS_corr')); % S
    IS = rescaleNPX(session, IS);
    
    %% Load Seizure Events
    try
        seizures = dlmread(strcat(path,'/',sessionID,'.seizures'));
        seizures = rescaleNPX(session, seizures);
    catch
        seizures = struct();
    end


    %% Load Slow-Wave Sleep (SWS) periods
    % Possibility to have an absence of a file or an empty file. in case of
    % a XOR empty file between sws and rem, we need to affect the default
    % value empty and not NaN.
    try
        swsIntervals=dlmread(strcat(path,'/',sessionID,'.sws'));
    catch
        swsIntervals=[];
    end
    try 
        remIntervals=dlmread(strcat(path,'/',sessionID,'.rem'));
    catch
        remIntervals=[];
    end
    sleepIntervals = [swsIntervals;remIntervals];

    if verbose 
        if isempty(sleepIntervals)
            fprintf("No sleep for %s",char(strcat(path,'/',sessionID)))
        elseif isempty(swsIntervals)
            fprintf(".sws file absent or empty for %s",char(strcat(path,'/',sessionID)))

        elseif isempty(remIntervals)
            fprintf(".rem file absent or empty for %s",char(strcat(path,'/',sessionID)))
        end
    end        

    if ~isempty(sleepIntervals)
        sleepIntervals = sortrows(sleepIntervals, 1);
        sleep = rescaleNPX(session, sleepIntervals);

        if ~isempty(swsIntervals)
            sws = rescaleNPX(session, swsIntervals);
        end

        % Erase micro awakness if needed
        if tsMicroAwake > 0 % Only if valid and not null treshold
            sleep.data = cleaned_intervals(sleep.data, tsMicroAwake);  % Convert micro-awakenings (<tsMicroAwake in s) into continuous sleep periods.
        end
        
        % Debugging
        checkIntervals(sessionID, sleep.data, 60*60+30*60,20) % threshold to output warning (1h30 for sleep, 9000s=2h30 for wake below)

        % WAKE INTERVALS - Thanks the sleep intervals and it's unrecorded intervals compute wake intervals
        windowTimestamps = [0 seconds(sleep.refStartAndStop(2)-sleep.refStartAndStop(1))];
        unrecInts = sleep.unrecInt;
        stateVect = sortrows([sleep.data;unrecInts]);
        wake.data = zeros(size(stateVect,1)+1,2);
        wake.data(1,1) = windowTimestamps(1);
        wake.data(end,2) = windowTimestamps(2);
        wake.data(2:end,1) = stateVect(:,2);
        wake.data(1:end-1,2) = stateVect(:,1);
        wake.refStartAndStop = sleep.refStartAndStop;
        wake.unrecInt = sleep.unrecInt;

         % Debugging
        checkIntervals(sessionID, wake.data, 2*60*60+30*60,20)
    else
        sleep = struct();
        wake = struct();
        sws = struct();
    end

        % Keep a track on sws intervals -> deduce rem with sws and sleep
        % (if tsMicroAwake > 0 then cleaned_intervals add some sleep time
        % due to "false awaking" -> labelled as rem because LFP of rem closer to awaken)


    %% Extra Features -> TBD by Sara
    % is_type={dlmread(strcat(path,'/',sessionID,'.IS_type'))}; % timestamps in SECONDS (divide by 3600 if you want HOURS as in IS)
    % is_type= {NPXrescaleZT_v3(session, is_type, ZT)};
    ISType= struct();


    % waveforms={dlmread(strcat(path,'/',sessionID,'.cortex_IS_waveforms_corr'))};
    ISWaveforms = struct();
    


    % try
    %     dRipples=dlmread(strcat(path,'/',sessionID,'.dRipples_corr'));
    %     dRipples={NPXrescaleZT_v2(session, dRipples(:,[1,3]), 0)};
    % catch
    %     dRipples=NaN;
    % end

    %% ZeitGeber rescaling -> to complete with extra features if needed
    if ~isnan(ZT)
        IS = rescaleZT(IS,ZT);
        if ~isempty(fieldnames(sleep))
            sleep = rescaleZT(sleep,ZT);
            wake = rescaleZT(sleep,ZT);
        end
        if ~isempty(fieldnames(sws))
            sws = rescaleZT(sws,ZT);
        end
        if ~isempty(fieldnames(seizures))
            seizures = rescaleZT(seizures,ZT);
        end
    else  
        IS.ZT = NaN;
        if ~isempty(fieldnames(sleep)) 
            sleep.ZT = NaN;
            wake.ZT = NaN;
        end
        if ~isempty(fieldnames(sws))
            sws.ZT = NaN;
        end
        if ~isempty(fieldnames(seizures))
            seizures.ZT = NaN;
        end
    end
    
end



