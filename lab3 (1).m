function result = lab3(varargin)
% Lab 3. Full and reduced Mamdani rule bases.
% Run LAB3 for the plots and the R2014 Rule Editor.
% LAB3(false) runs the calculations and plots without opening the editor.
if nargin == 2
    % The legacy Fuzzy Logic Toolbox calls this same file as lab3(x,params)
    % for a custom membership function. No extra M-file is needed.
    result = interpolateTerm(varargin{1}, varargin{2});
    return
end
openEditor = nargin == 0 || varargin{1};

% Parameters are the eight membership values at the interval centers
% of each variable in Lab 2. Rows of a block are: Low, Medium and High.
ranges = [60 100; 2 10; 5 9; 0 100];
nodes = {
    [1 0.7135254915624212 0.7135254915624212 0.410218024111047 0.15450849718747373 0 0 0;
          0 0.21654236465910048 0.21654236465910048 0.8106427859674035 1 0.5749171248427078 0 0;
          0 0 0 0 0 0.5784142923389222 1 1]
    [1 0.7135254915624212 0.410218024111047 0 0 0 0 0;
          0 0.15450849718747373 0.5784142923389222 1 0.5784142923389222 0.410218024111047 0 0;
          0 0 0 0 0.410218024111047 0.5784142923389222 1 1]
    [1 0.7135254915624212 0.7135254915624212 0 0 0 0 0;
          0 0.15450849718747373 0.15450849718747373 1 1 0.5784142923389222 0 0;
          0 0 0 0 0 0.410218024111047 1 1]
    [1 1 0.7135254915624212 0.11435677693908453 0 0 0 0;
          0 0 0.15450849718747373 0.8187293044088437 1 0.7135254915624212 0 0;
          0 0 0 0 0 0.15450849718747373 1 1]
};
params = cell(4, 3);
for k = 1:4
    for t = 1:3
        params{k, t} = [ranges(k, :), nodes{k}(t, :)];
    end
end

% All 27 combinations of input terms. A 0 in a reduced rule omits cargo.
fullRules = zeros(27, 6);
reducedRules = zeros(27, 6);
n = 0;
m = 0;
for axle = 1:3
    for tyre = 1:3
        cargo = (1:3).';
        group = [cargo, axle*ones(3, 1), tyre*ones(3, 1)];

        % Classroom rule policy: high axle load targets the high pressure
        % term. Below-target pressure has the larger penalty.
        target = 2 + (axle == 3);
        risk = 0.20*(cargo-1) + 0.40*(axle-1) + ...
               0.90*max(target-tyre, 0) + 0.45*max(tyre-target, 0);
        outcome = 3 - min(2, round(risk));
        group = [group, outcome, ones(3, 2)];
        fullRules(n+1:n+3, :) = group;
        n = n + 3;
        if all(outcome == outcome(1))
            m = m + 1;
            reducedRules(m, :) = [0 axle tyre outcome(1) 1 1];
        else
            reducedRules(m+1:m+3, :) = group;
            m = m + 3;
        end
    end
end
reducedRules = reducedRules(1:m, :);
for k = 1:n
    exact = bsxfun(@eq, reducedRules(:, 1:3), fullRules(k, 1:3));
    match = all(exact | reducedRules(:, 1:3) == 0, 2);
    assert(nnz(match) == 1 && reducedRules(match, 4) == fullRules(k, 4), ...
           'A term combination is missing or has conflicting conclusions.');
end
fprintf('Rules: %d full, %d reduced. All 27 term combinations covered.\n', n, m);

% Product AND and implication, sum aggregation, sampled centroid.
y = (0:100).';
outMu = zeros(numel(y), 3);
for k = 1:3
    outMu(:, k) = interpolateTerm(y, params{4, k});
end
areas = sum(outMu, 1).';
moments = (y.' * outMu).';
compute = @(x, r) evaluateRules(x, r, params, areas, moments);

% Cases 1-3: EPA track data (PMC6198682, Table 1 and axle-load text).
% Drive and trailer totals are split evenly across their two axles.
% Cases 4-5: Dunlop LI 160 single-tyre axle capacity table.
% Cargo 80 kg is held constant for the test. Actual safety is not measured.
tests = [80 5.37 7.6;
          80 6.07 6.9;
          80 6.145 6.9;
          80 8.3 8;
          80 8.3 8.5];
capacity = [NaN; NaN; NaN; 8.19; 8.60];
fullTest = compute(tests, fullRules);
shortTest = compute(tests, reducedRules);
result = [tests, capacity, fullTest, shortTest, shortTest-fullTest];
fprintf('\nCargo kg  Axle t  Pressure bar  Capacity t     Full  Reduced   Change\n');
fprintf('%8.0f %8.3f %13.2f %11.3f %8.2f %8.2f %8.2f\n', result.');
assert(fullTest(5) > fullTest(4) && shortTest(5) > shortTest(4), ...
       'The catalog comparison has the wrong direction.');

% Compare a dense surface slice using common axes for both bases.
[cargo, axle] = meshgrid(linspace(60, 100, 41), linspace(2, 10, 41));
pressure = 7;
input = [cargo(:), axle(:), pressure*ones(numel(cargo), 1)];
zFull = reshape(compute(input, fullRules), size(cargo));
zReduced = reshape(compute(input, reducedRules), size(cargo));
delta = zReduced - zFull;
fprintf('Pressure slice: %.1f bar. Max. difference: %.3f points.\n', ...
        pressure, max(abs(delta(:))));
fig = figure('Name', 'Lab 3 - Rule surfaces', 'NumberTitle', 'off', ...
             'Color', 'w', 'Position', [60 90 1200 440]);
subplot(1, 3, 1);
surf(cargo, axle, zFull, 'EdgeColor', 'none');
labelSurface(sprintf('Full base: %d rules', n), pressure);
subplot(1, 3, 2);
mesh(cargo, axle, zReduced);
labelSurface(sprintf('Reduced base: %d rules', m), pressure);
subplot(1, 3, 3);
contourf(cargo, axle, delta, 16, 'LineColor', 'none');
title('Reduced minus full');
xlabel('Cargo weight (kg)');
ylabel('Axle load (t)');
colorbar;

folder = fullfile(fileparts(mfilename('fullpath')), 'graphs');
if ~exist(folder, 'dir')
    mkdir(folder);
end
print(fig, fullfile(folder, 'lab3_rule_surfaces.png'), '-dpng', '-r150');

% The positional NEWFIS syntax is the form used by MATLAB R2014.
if exist('newfis', 'file') && license('test', 'Fuzzy_Toolbox')
    fisFull = makeFIS('TruckSafetyFull', fullRules, ranges, params);
    fisReduced = makeFIS('TruckSafetyReduced', reducedRules, ranges, params);
    writefis(fisFull, fullfile(folder, 'TruckSafetyFull'));
    writefis(fisReduced, fullfile(folder, 'TruckSafetyReduced'));
    % The 2014 form of EVALFIS takes inputs first and the FIS second.
    try
        nativeFull = evalfis(tests, fisFull, numel(y));
        nativeReduced = evalfis(tests, fisReduced, numel(y));
        difference = max(abs([nativeFull(:) - fullTest; nativeReduced(:) - shortTest]));
        fprintf('Maximum difference versus Toolbox EVALFIS: %.3g.\n', difference);
    catch err
        warning('Could not compare against Toolbox EVALFIS: %s', err.message);
    end
    if openEditor
        ruleedit(fisFull);
    end
else
    fprintf('Fuzzy Logic Toolbox is needed to open the R2014 Rule Editor.\n');
end
end

function mu = interpolateTerm(x, p)
% The FIS stores just 10 parameters: the two bounds and eight node values.
edges = linspace(p(1), p(2), 9);
centers = (edges(1:end-1) + edges(2:end)) / 2;
knots = [p(1), centers, p(2)];
levels = p(3:10);
levels = [levels(1), reshape(levels, 1, []), levels(end)];
x = min(p(2), max(p(1), x));
mu = interp1(knots, levels, x, 'linear');
mu = reshape(mu, size(x));
end

function y = evaluateRules(input, rules, params, areas, moments)
% Accumulate firing strengths by consequent. Product-sum aggregation makes
% the area and first moment linear in these three accumulated strengths.
n = size(input, 1);
terms = zeros(n, 3, 3);
for k = 1:3
    for t = 1:3
        terms(:, t, k) = interpolateTerm(input(:, k), params{k, t});
    end
end
totals = zeros(n, 3);
for row = 1:size(rules, 1)
    strength = rules(row, 5) * ones(n, 1);
    for k = 1:3
        t = rules(row, k);
        if t ~= 0
            strength = strength .* terms(:, t, k);
        end
    end
    out = rules(row, 4);
    totals(:, out) = totals(:, out) + strength;
end
denominator = totals * areas;
assert(all(denominator > 0), 'No rule is active at one or more points.');
y = (totals * moments) ./ denominator;
end

function labelSurface(caption, pressure)
title(sprintf('%s, %.1f bar', caption, pressure));
xlabel('Cargo weight (kg)');
ylabel('Axle load (t)');
zlabel('Safety score (%)');
axis([60 100 2 10 0 100]);
caxis([0 100]);
view(-42, 30);
grid on;
end

function fis = makeFIS(name, rules, ranges, params)
fis = newfis(name, 'mamdani', 'prod', 'max', 'prod', 'sum', 'centroid');
names = {'CargoWeight', 'AxleLoad', 'TyrePressure', 'Safety'};
terms = {'Low', 'Medium', 'High'};
for k = 1:4
    if k <= 3
        kind = 'input';
        index = k;
    else
        kind = 'output';
        index = 1;
    end
    fis = addvar(fis, kind, names{k}, ranges(k, :));
    for t = 1:3
        fis = addmf(fis, kind, index, terms{t}, 'lab3', params{k, t});
    end
end
fis = addrule(fis, rules);
end
