% Load the data
data = readmatrix('scaled_WT_points.csv');  % Make sure this CSV file is in your path
x = data(:,1);
y = data(:,2);
n = length(x);

% Define 100 box sizes in log10 scale
k = 0:99;
log_epsilons = -1 - 0.025 * k;
epsilons = 10.^log_epsilons;
N_eps = zeros(size(epsilons));

% Loop over each epsilon
for idx = 1:length(epsilons)
    epsilon = epsilons(idx);
    delta = epsilon / 2;
    Z = [];

    % Interpolate each segment
    for i = 1:n-1
        dx = x(i+1) - x(i);
        dy = y(i+1) - y(i);
        Li = sqrt(dx^2 + dy^2);
        Mi = max(ceil(Li / delta), 1);
        t = linspace(0, 1, Mi+1);
        xi = x(i) + t * dx;
        yi = y(i) + t * dy;
        Z = [Z; xi(:), yi(:)];
    end

    % Remove duplicate points
    Z = unique(Z, 'rows');

    % Compute midpoints and assign to boxes
    midpoints = (Z(1:end-1,:) + Z(2:end,:)) / 2;
    box_indices = floor(midpoints / epsilon);
    [~, unique_idx] = unique(box_indices, 'rows');
    N_eps(idx) = length(unique_idx);
end

% Plot log-log graph
figure;
loglog(epsilons, N_eps, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 10);
xlabel('log_{10}(\epsilon)');
ylabel('log_{10}(N(\epsilon))');
title('Dense path sampling wt');
grid on;

% Overlay log-log ticks
set(gca, 'XScale', 'log', 'YScale', 'log');
