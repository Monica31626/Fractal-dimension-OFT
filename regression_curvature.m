% Load data with headers using readmatrix
data = readmatrix('scaled_WT_points.csv');  % [x, y] columns
x = data(:, 1);
y = data(:, 2);
n = length(x);

% Parameters
S_base = 1;
alpha = 100;

% Step 1: Curvature-based oversampling
Z = []; % Oversampled points
for i = 2:n-1
    xi_prev = [x(i-1), y(i-1)];
    xi = [x(i), y(i)];
    xi_next = [x(i+1), y(i+1)];

    % Compute signed-area curvature
    num = 2 * abs(det([xi_next - xi; xi - xi_prev]));
    denom = norm(xi_next - xi) * norm(xi - xi_prev) * norm(xi_next - xi_prev);

    if denom == 0
        kappa = 0;
    else
        kappa = num / denom;
    end

    Si = S_base + floor(alpha * kappa);
    for s = 1:Si
        t = s / (Si + 1);
        gi = xi + t * (xi_next - xi);
        Z = [Z; gi];
    end
end
Z = [Z; x(end), y(end)];

% Step 2: Box counting
epsilons = 10.^(-0.025 * (0:400));
N_eps = zeros(size(epsilons));
for i = 1:length(epsilons)
    eps = epsilons(i);
    m = floor(Z / eps);
    unique_boxes = unique(m, 'rows');
    N_eps(i) = size(unique_boxes, 1);
end

% Step 3: Filter scales based on regression line
log_eps = log10(epsilons);
log_N = log10(N_eps);

N_max = max(N_eps);
max_idx = find(N_eps == N_max, 1);  % first occurrence
reg_range = 1:max_idx;
p = polyfit(-log_eps(reg_range), log_N(reg_range), 1);
a = p(1); b = p(2);

mask = log_N > (-a * log_eps + b);
S1_idx_all = find(mask);

% Find longest consecutive subset of indices
if length(S1_idx_all) < 3
    error('Too few scales passed the regression line. Check parameters or data.');
end

% Detect breaks in consecutiveness
d = diff(S1_idx_all);
if isempty(d)
    S1_idx = S1_idx_all;
else
    ends = [0; find(d > 1); length(S1_idx_all)];
    max_len = 0; best_start = 1;

    for j = 1:length(ends)-1
        len = ends(j+1) - ends(j);
        if len > max_len
            max_len = len;
            best_start = ends(j) + 1;
        end
    end

    S1_idx = S1_idx_all(best_start:best_start + max_len - 1);
end

% Step 4: Apply TLS on each segment and compute SSE
S1_log_eps = log_eps(S1_idx);
S1_log_N = log_N(S1_idx);
K = length(S1_idx);
SSE = zeros(1, K - 3);

for i = 2:K-2
    x0 = -S1_log_eps(1:i);
    y0 = S1_log_N(1:i);
    x1 = -S1_log_eps(i+1:end);
    y1 = S1_log_N(i+1:end);

    p0 = polyfit(x0, y0, 1);
    p1 = polyfit(x1, y1, 1);
    err0 = sum((y0 - polyval(p0, x0)).^2);
    err1 = sum((y1 - polyval(p1, x1)).^2);
    SSE(i) = err0 + err1;
end

% Find optimal elbow
[~, i_star] = min(SSE);
elbow_eps = epsilons(S1_idx(i_star));

% Step 5: Plot
figure;
plot(-log_eps(S1_idx(2:K-2)), SSE(2:K-2), 'b.-', 'LineWidth', 1.5); hold on;
xline(-log10(elbow_eps), 'r--', 'LineWidth', 1.5);
xlabel('-log_{10}(\epsilon)');
ylabel('Total SSE');
title('Regression-Based TLS Elbow Detection (Curvature FD)');
legend('SSE', 'Elbow Point');
grid on;
