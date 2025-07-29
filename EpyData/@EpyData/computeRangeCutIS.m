function this = computeRangeCutIS(this)
%COMPUTERANGECUTIS Sets manual time cutoffs for IS count to avoid edge artifacts
%
% This method sets the field `rangeCutIS` of the EpyData object to define
% per-rat time limits (in seconds) up to which interictal spike (IS) counts 
% are valid and free from edge effects.
%
% The cutoff values were identified graphically using the function
% `plotPSDPerRat` by observing IS count artifacts at the end of some sessions.
% For rats with no visible edge effect, the value is left empty ([]).
%
% The function distinguishes between Neuropixels and Telemetry recordings
% based on the sampling frequency, and assigns the appropriate cutoff values.

arguments
    this (1,1) EpyData
end

%% Create a structure storing cutoff times per rat (in seconds)

this.rangeCutIS = zeros(numel(this.ratID),1);
timestampsRangeStruct = struct();

% -------- Neuropixels (NPX) --------
% Rats for which no end-of-session artifacts were detected → no cutoff needed
npx_ids = [596, 597, 671, 673, 674, 794, 796];
for id = npx_ids
    timestampsRangeStruct.NPX.(sprintf('rat%d', id)) = [];
end

% -------- Telemetry (TM) --------
% Some rats showed edge artifacts → set cutoff times (in hours)
tm_data = {
    677,   79.7917;       
    982,   80.375;
    986,   47;    
    1065,  76.375;
    1069,  77.7083;
    1070,  74;       % Not taken in account : eeg 37, 38, 39, 40, 41, 42, 43, 44
    1074,  73;         % 74.9586 days (former)
    1212,  [];        % No artifact
    1213,  54.8333;
    1214,  54.875;
    1218,  54.875;
    1220,  28.875;
    1224,  62.875;
};

% Convert times to seconds and store in struct
for k = 1:size(tm_data,1)
    id = tm_data{k,1};
    val = tm_data{k,2};
    timestampsRangeStruct.TM.(sprintf('rat%d', id)) = val * 24 * 3600; % hours → seconds
end

%% Identify recording modality based on sampling frequency
if this.samplingFreq == 30000
    modality = "NPX"; % Neuropixels
elseif this.samplingFreq == 512
    modality = "TM";  % Telemetry
else
    error('Unsupported sampling frequency: %d', this.samplingFreq);
end

%% Assign the cutoff values to the object for each rat
% 'ratID' should match the 'rat###' format used in the struct fields
for i = 1:numel(this.ratID)
    if ~isfield(timestampsRangeStruct.(modality), sprintf('rat%d', this.ratID(i))) || isempty(timestampsRangeStruct.(modality).(sprintf('rat%d', this.ratID(i)))) 
        this.rangeCutIS(i,1) = seconds(this.refDatesTimestamps(i,2)-this.refDatesTimestamps(i,1));
    else 
        this.rangeCutIS(i,1) = timestampsRangeStruct.(modality).(sprintf('rat%d', this.ratID(i)));
    end
end

end
