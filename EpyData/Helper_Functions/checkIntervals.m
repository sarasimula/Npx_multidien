function checkIntervals(sessionID, intervals, durationThreshold, countThreshold)
% checkIntervals Vérifie les intervalles selon des critères de durée et de nombre.
%
% Syntaxe :
%   checkIntervals(sessionID, intervals, durationThreshold, countThreshold)
%
% Entrées :
%   sessionID         : string ou char, identifiant de la session
%   intervals         : matrice Nx2, avec [début, fin] des intervalles
%   durationThreshold : durée minimale acceptable d’un intervalle (en secondes)
%   countThreshold    : nombre minimal d’intervalles requis
%
% Affiche un warning si certains intervalles sont trop longs ou s'il y en a trop peu.

    % Sanity check
    if size(intervals,2) ~= 2
        error('La matrice d''intervalles doit avoir deux colonnes [start, end].');
    end

    % Durées des intervalles
    durations = intervals(:,2) - intervals(:,1);

    % Vérifie les durées trop longues
    idx_long = find(durations >= durationThreshold);
    if ~isempty(idx_long)
        dur_long = durations(idx_long);
        warning('Session %s : %d intervalle(s) ont une durée >= %g s. (Indices : %s with %s s)', ...
            sessionID, numel(idx_long), durationThreshold, mat2str(idx_long'), mat2str(dur_long'));
    end

    % Vérifie les durées négatives
    idx_long = find(durations < 0);
    if ~isempty(idx_long)
        warning('Session %s : %d intervalle(s) ont une durée < 0 s. (Indices : %s)', ...
            sessionID, numel(idx_long), mat2str(idx_long'));
    end

    % Vérifie le nombre d’intervalles
    if size(intervals,1) <= countThreshold
        warning('Session %s : not enough sleep intervals (seulement %d, seuil = %d).', ...
            sessionID, size(intervals,1), countThreshold);
    end
end
