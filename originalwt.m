% Load the trajectory data (assumed to be 2D)
data = readmatrix('scaled_WT_points.csv');
X = data;  % X(:,1) = x-coordinates, X(:,2) = y-coordinates
n = size(X, 1);

% Define box sizes: log10(epsilon_k) = -1 - 0.025*k, for k = 0 to 99
k_vals = 0:99;
log_eps = -1 - 0.025 * k_vals;
epsilons = 10.^log_eps;

% Initialize N(epsilon) vector
N_eps = zeros(length(epsilons), 1);

% Loop over each epsilon
for i = 1:length(epsilons)
    epsilon = epsilons(i);

    % Compute the index of the ε-box each point falls into
    box_indices = floor(X / epsilon);  % (n x 2) array of box coordinates

    % Count the number of unique boxes
    unique_boxes = unique(box_indices, 'rows');
    N_eps(i) = size(unique_boxes, 1);
end

% Log-log values
log_N = log10(N_eps);
log_epsilons = log10(epsilons);

% Linear regression for FD estimate
fit_range = 20:80;  % Choose scale range that looks linear
p = polyfit(log_epsilons(fit_range), log_N(fit_range), 1);
DB = -p(1);

% Display result
fprintf('Estimated Box-Counting Dimension D_B = %.4f\n', DB);

% Plot log-log graph
figure;
plot(log_epsilons, log_N, 'bo-', 'LineWidth', 1.3, 'MarkerSize', 6); hold on;
plot(log_epsilons(fit_range), polyval(p, log_epsilons(fit_range)), 'r--', 'LineWidth', 1.5);
xlabel('log_{10}(\epsilon)');
ylabel('log_{10}(N(\epsilon))');
title('Original Box-Counting - WT');
legend('Data', sprintf('Linear Fit'), 'Location', 'SouthWest');
grid on;
