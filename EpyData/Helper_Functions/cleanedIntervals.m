function I_clean = cleanedIntervals(I, t)
% cleanedIntervals - Convert some micro state intervals (the exterior state of the matrix I) into the opposite state (the inner state of the matrix I).
%
% USAGE:
%   I_clean = cleanedIntervals(I, t)
%
% INPUTS:
%   I : Nx2 double, indicating intervals of a state. Colum 1 : beggining of
%   intervals. Column 2 : ending of intervals.
%       
%   t : 1x1 positive double, treshold below which micro-state intervals are merged into the other state intervals.
%
% OUTPUT:
%   I_clean : Mx2 double, new I vector with the conversion operated.

% Example :     
% I : Nx2 tableau des intervalles de sommeil par exemple
% t : seuil minimal de durée pour l'éveil à respecter (on veut enlever
% les micro éveils)

arguments
    I (:,2) double % string or char
    t (1,1) double {mustBeNonnegative}
end


% Initialisation
I_clean = I(1,:); % First interval

for k = 2:size(I,1)
    % Duration of the "exterior" state between I(k-1) and I(k)
    sleep_duration = I(k,1) - I_clean(end,2);

    if sleep_duration < t % If strictly below 
        % Micro-state : merge the two inners state intervals
        I_clean(end,2) = I(k,2);
    else
        % Duration of exterior state enough : new segment
        I_clean = [I_clean; I(k,:)];
    end
end

end