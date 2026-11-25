function lab2
% Lab 2. Four membership plots from two rounds of pairwise comparisons.
% All functions and data used by this file are included below.
names = {'Cargo weight', 'Axle load', 'Tyre pressure', 'Driving safety'};
units = {'Cargo weight (kg)', 'Axle load (t)', 'Tyre pressure (bar)', ...
         'Driving safety (%)'};
edges = {60:5:100, 2:10, 5:0.5:9, 0:12.5:100};
files = {'lab2_cargoweight.png', 'lab2_axleload.png', ...
         'lab2_tirepressure.png', 'lab2_safety.png'};
terms = {'Low', 'Medium', 'High'};
colors = [0.27 0.45 0.77; 0.93 0.49 0.19; 0.44 0.68 0.28];
styles = {'-', '--', '-.'};

% Each block contains 5 experts x 3 terms for 8 intervals.
votes = {
    [1 1 1 0 0 0 0 0;
          0 0 0 1 1 0 0 0;
          0 0 0 0 0 1 1 1;
          1 1 1 0 0 0 0 0;
          0 0 0 1 1 0 0 0;
          0 0 0 0 0 1 1 1;
          1 1 1 1 0 0 0 0;
          0 0 0 0 1 1 0 0;
          0 0 0 0 0 0 1 1;
          1 1 1 1 1 0 0 0;
          0 0 0 0 0 1 0 0;
          0 0 0 0 0 0 1 1;
          1 0 0 0 0 0 0 0;
          0 1 1 1 1 0 0 0;
          0 0 0 0 0 1 1 1]
    [1 1 0 0 0 0 0 0;
          0 0 1 1 1 0 0 0;
          0 0 0 0 0 1 1 1;
          1 1 0 0 0 0 0 0;
          0 0 1 1 0 0 0 0;
          0 0 0 0 1 1 1 1;
          1 1 1 0 0 0 0 0;
          0 0 0 1 1 1 0 0;
          0 0 0 0 0 0 1 1;
          1 1 1 0 0 0 0 0;
          0 0 0 1 1 1 0 0;
          0 0 0 0 0 0 1 1;
          1 0 0 0 0 0 0 0;
          0 1 1 1 0 0 0 0;
          0 0 0 0 1 1 1 1]
    [1 1 1 0 0 0 0 0;
          0 0 0 1 1 0 0 0;
          0 0 0 0 0 1 1 1;
          1 1 1 0 0 0 0 0;
          0 0 0 1 1 1 0 0;
          0 0 0 0 0 0 1 1;
          1 1 1 0 0 0 0 0;
          0 0 0 1 1 1 0 0;
          0 0 0 0 0 0 1 1;
          1 1 1 0 0 0 0 0;
          0 0 0 1 1 1 0 0;
          0 0 0 0 0 0 1 1;
          1 0 0 0 0 0 0 0;
          0 1 1 1 1 0 0 0;
          0 0 0 0 0 1 1 1]
    [1 1 1 1 0 0 0 0;
          0 0 0 0 1 1 0 0;
          0 0 0 0 0 0 1 1;
          1 1 1 0 0 0 0 0;
          0 0 0 1 1 1 0 0;
          0 0 0 0 0 0 1 1;
          1 1 1 0 0 0 0 0;
          0 0 0 1 1 1 0 0;
          0 0 0 0 0 0 1 1;
          1 1 1 0 0 0 0 0;
          0 0 0 1 1 1 0 0;
          0 0 0 0 0 0 1 1;
          1 1 0 0 0 0 0 0;
          0 0 1 1 1 0 0 0;
          0 0 0 0 0 1 1 1]
};

folder = fullfile(fileparts(mfilename('fullpath')), 'graphs');
if ~exist(folder, 'dir')
    mkdir(folder);
end
for variable = 1:4
    bounds = edges{variable};
    centers = (bounds(1:end-1) + bounds(2:end)) / 2;
    mu = comparePanels(votes{variable});
    fprintf('\n%s\n', names{variable});
    fprintf('        Center        Low     Medium       High\n');
    fprintf('%14.2f %10.3f %10.3f %10.3f\n', [centers.', mu].');

    % Extend the end values to the limits, then interpolate linearly.
    knots = [bounds(1), centers, bounds(end)];
    x = linspace(bounds(1), bounds(end), 401).';
    curves = interp1(knots, mu([1 1:8 8], :), x, 'linear');
    fig = figure('Name', ['Lab 2 - ' names{variable}], ...
                 'NumberTitle', 'off', 'Color', 'w');
    lines = plot(x, curves, 'LineWidth', 1.8);
    for term = 1:3
        set(lines(term), 'Color', colors(term, :), 'LineStyle', styles{term});
    end
    title(names{variable});
    xlabel(units{variable});
    ylabel('Membership degree');
    set(gca, 'XTick', bounds);
    axis([bounds(1) bounds(end) 0 1.05]);
    legend(terms, 'Location', 'SouthOutside', 'Orientation', 'horizontal');
    grid on;
    print(fig, fullfile(folder, files{variable}), '-dpng', '-r150');
end
end

function mu = comparePanels(votes)
% Experts 1-3 and 4-5 form the two panels. Convert support p to q = 1 + 8p.
mu = zeros(8, 3);
for term = 1:3
    t = votes(term:3:end, :);
    q1 = 1 + 8 * mean(t(1:3, :), 1).';
    q2 = 1 + 8 * mean(t(4:5, :), 1).';
    A = bsxfun(@rdivide, q1, q1.');
    B = bsxfun(@rdivide, q2, q2.');
    joint = sqrt(A .* B);
    weights = joint(:, 1) / sum(joint(:, 1));

    % Remove the baseline for zero support and scale the peak to one.
    baseline = 1 / sum(sqrt(q1 .* q2));
    if any(t(:))
        mu(:, term) = max(0, weights - baseline) / (max(weights) - baseline);
    end
end
end
