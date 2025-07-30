% Load CSV data (assuming two columns: x and y)
data = readmatrix('wt_stress_points.xlsx');
n = size(data, 1);

% Define epsilon values: log10(epsilon_k) = -1 - 0.025*k, for k = 0,...,99
k_vals = 0:99;
log_epsilons = -1 - 0.025 * k_vals;
epsilons = 10.^log_epsilons;

% Define delta values: δ = 1,...,100
deltas = 1:100;

% Initialize N(epsilon, delta) matrix
N_eps_delta = zeros(length(epsilons), length(deltas));

% Loop over delta values
for j = 1:length(deltas)
    delta = deltas(j);
    m_delta = floor((n - 1) / delta);

    % Get resampled points zk = x_{k*delta + 1}
    indices = delta * (0:m_delta) + 1;
    indices(indices > n) = []; % Ensure valid indices
    zk = data(indices, :);

    % Loop over epsilon values
    for i = 1:length(epsilons)
        epsilon = epsilons(i);
        rk_sum = 0;

        for k = 1:(size(zk, 1) - 1)
            dist = norm(zk(k+1,:) - zk(k,:));
            rk = ceil(dist / epsilon);
            rk_sum = rk_sum + rk;
        end

        N_eps_delta(i, j) = rk_sum;
    end
end

% ----- Plotting and FD estimation for delta = 10 (index 10) -----
delta_idx = 10; % delta = 10
logN = log10(N_eps_delta(:, delta_idx));
p = polyfit(log_epsilons', logN, 1);
FD_estimate = -p(1); % Slope is -FD

% Plotting
figure;
plot(log_epsilons, logN, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 10);
hold on;
plot(log_epsilons, polyval(p, log_epsilons), 'r--', 'LineWidth', 1.5);
xlabel('log_{10}(\epsilon)');
ylabel('log_{10}(N(\epsilon, \delta))');
title(['Temporal-spatial Divider-wt+stress']);
legend('Data', 'Linear fit', 'Location', 'SouthWest');
grid on;
