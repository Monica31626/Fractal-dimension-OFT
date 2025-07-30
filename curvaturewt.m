% Load scaled data
data = readmatrix('scaled_WT_points.csv'); % Ensure this file is in your working directory
x = data(:,1);
y = data(:,2);
n = length(x);

% Parameters
Sbase = 1;       % Minimum number of interpolated points
alpha = 10;      % Curvature sensitivity
Z = [];          % Oversampled trajectory

for i = 2:n-1
    % Define three consecutive points
    xi_prev = [x(i-1), y(i-1)];
    xi = [x(i), y(i)];
    xi_next = [x(i+1), y(i+1)];

    % Compute curvature kappa_i
    v1 = xi_next - xi;
    v2 = xi - xi_prev;
    cross = v1(1)*v2(2) - v1(2)*v2(1);  % 2D cross product magnitude
    norm1 = norm(v1);
    norm2 = norm(v2);
    norm3 = norm(xi_next - xi_prev);
    kappa = 2 * abs(cross) / (norm1 * norm2 * norm3 + eps); % add eps to avoid divide-by-zero

    % Determine oversampling count
    Si = Sbase + floor(alpha * kappa);

    % Interpolate between xi and xi+1
    segment = [xi]; % start with xi
    for s = 1:Si
        t = s / (Si + 1);
        interp_point = xi + t * (xi_next - xi);
        segment = [segment; interp_point];
    end
    Z = [Z; segment];  % accumulate
end

% Append final point
Z = [Z; x(end), y(end)];

% Box-counting algorithm
epsilons = logspace(-3.5, -1, 100); % 100 log-spaced box sizes
N_eps = zeros(size(epsilons));

for j = 1:length(epsilons)
    eps = epsilons(j);
    % Discretize space into grid
    grid_x = floor(Z(:,1) / eps);
    grid_y = floor(Z(:,2) / eps);
    boxes = unique([grid_x, grid_y], 'rows');
    N_eps(j) = size(boxes,1);
end

% Estimate fractal dimension
log_eps = log10(epsilons);
log_N = log10(N_eps);
p = polyfit(log_eps, log_N, 1);
D_B = -p(1);

% Display result
fprintf('Estimated Fractal Dimension D_B = %.4f\n', D_B);

% Optional: Plot
figure;
plot(log_eps, log_N, 'bo-'); hold on;
plot(log_eps, polyval(p, log_eps), 'r--');
xlabel('log_{10}(\epsilon)');
ylabel('log_{10}(N(\epsilon))');
title(sprintf('curvature based Estimate wt'));
grid on;
