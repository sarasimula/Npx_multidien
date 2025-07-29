function h = plotComparativeViolin(data, xlabels)
% plotComparativeViolin - Compare two group distributions with violin plot, 
% box plots and scatters (+ stats).
%
% Inputs:
%   data    : (2,1) cell array of numerical vectors (group 1 and 2)
%   xlabels : (1,2) cell array of strings, x-axis labels for each group
%
% Output:
%   h : handle to the violin plot
%
% Uses Daviolinplot function form toolbox

% Violin plot
h = daviolinplot(data, ...
    'groups', [ones(size(data{1})); 2*ones(size(data{2}))], ...
    'xtlabels', xlabels, ...
    'colors', lines(2), ...
    'box', 2, 'scatter', 1, ...
    'jitter', 1, 'violin', 'half');

% Normality test (Lilliefors)
normal1 = ~lillietest(data{1});
normal2 = ~lillietest(data{2});

% Choose statistical test and compute effect size
if normal1 && normal2
    % If both distributions are normal, use a t-test and Cohen's d.
    % Cohen's d measures standardized mean difference, assuming equal variances.
    [hyp, p, ~, stats] = ttest2(data{1}, data{2});
    test_name = 't-test';
    effect_size = computeCohensD(data{1}, data{2});
    effect_label = sprintf('Cohen''s d = %.2f', effect_size);
else
    % If at least one is non-normal, use ranksum and Cliff’s delta.
    % Cliff’s δ quantifies the probability that a value in one group
    % is larger than in the other, minus the reverse.
    [p, hyp, stats] = ranksum(data{1}, data{2});
    test_name = 'ranksum';
    effect_size = computeCliffsDelta(data{1}, data{2});
    effect_label = sprintf('Cliff''s δ = %.2f', effect_size);
end

% Annotate the plot with statistical result
hold on;
y_max = max([data{1}; data{2}]);
y_annot = y_max + 0.05 * range([data{1}; data{2}]);
plot([1 2], [y_annot y_annot], 'k', 'LineWidth', 1.5);
text(1.5, y_annot + 0.1 * range([data{1}; data{2}]), ...
    sprintf('%s, h = %d, p = %.3g, %s', test_name, hyp, p, effect_label), ...
    'HorizontalAlignment', 'center', 'FontSize', 10);
end

function d = computeCohensD(x1, x2)
% computeCohensD - Compute Cohen's d for two independent samples.
% -> Cohen's d assumes normally distributed, continuous data with
%    similar variance and interprets difference in terms of standard deviations.
%
% Inputs:
%   x1, x2 : numerical vectors (group 1 and 2)
%
% Output:
%   d : effect size (Cohen’s d)

n1 = numel(x1); n2 = numel(x2);
s1 = var(x1); s2 = var(x2);
s_pooled = sqrt(((n1 - 1)*s1 + (n2 - 1)*s2) / (n1 + n2 - 2));
d = (mean(x1) - mean(x2)) / s_pooled;
end

function delta = computeCliffsDelta(x1, x2)
% computeCliffsDelta - Compute Cliff's delta (non-parametric effect size).
% -> Cliff’s delta is a non-parametric effect size measuring dominance
%    (probability one group is greater than another), suitable for skewed or ordinal data.
%
% Inputs:
%   x1, x2 : numerical vectors (group 1 and 2)
%
% Output:
%   delta : Cliff’s δ, ranging from -1 to 1

n1 = numel(x1);
n2 = numel(x2);
count = 0;

for i = 1:n1
    count = count + sum(x1(i) > x2) - sum(x1(i) < x2);
end

delta = count / (n1 * n2);
end
