function run_dualTLS_original_boxcounting()
    % Step 1: Load data
    data = readtable('scaled_WT_points.csv');
    points = [data.x_prime, data.y_prime];  % Assumes columns: x_prime, y_prime

    % Step 2: Define epsilon values
    k = 0:800;
    epsilons = 10.^(-0.025 * k);
    num_eps = length(epsilons);
    
    % Step 3: Compute N(epsilon) for original box-counting algorithm
    N_eps = zeros(1, num_eps);
    for i = 1:num_eps
        eps = epsilons(i);
        box_indices = floor(points / eps);
        [~, unique_rows] = unique(box_indices, 'rows');
        N_eps(i) = length(unique_rows);
    end

    % Step 4: Log-log transformation
    log_eps_all = -log10(epsilons);           % x-axis
    log_N_all = log10(N_eps);                 % y-axis
    
    % Step 5: Plot Dual TLS with different gamma thresholds
    dualTLS_elbow(log_eps_all, log_N_all, N_eps);
end

function dualTLS_elbow(log_eps_all, log_N_all, N_eps)
    gamma_values = [0.01, 0.05];
    N_max = max(N_eps);
    figure; hold on;
    colors = ['b', 'r'];

    for idx = 1:length(gamma_values)
        gamma = gamma_values(idx);

        % Threshold filtering
        mask = (N_eps >= gamma * N_max) & (N_eps <= (1 - gamma) * N_max);
        log_eps = log_eps_all(mask);
        log_N = log_N_all(mask);
        K = length(log_eps);

        elbows = 3:K-2;
        sse = zeros(1, length(elbows));
        
        for j = 1:length(elbows)
            i = elbows(j);
            x0 = log_eps(1:i); y0 = log_N(1:i);
            x1 = log_eps(i+1:end); y1 = log_N(i+1:end);

            [m0, b0] = tls(x0, y0);
            [m1, b1] = tls(x1, y1);
            pred0 = m0 * x0 + b0;
            pred1 = m1 * x1 + b1;

            sse(j) = sum((y0 - pred0).^2) + sum((y1 - pred1).^2);
        end

        [~, min_idx] = min(sse);
        elbow_x = log_eps(elbows(min_idx));

        plot(log_eps(elbows), sse, colors(idx), 'LineWidth', 1.5); 
        xline(elbow_x, '--', ['γ = ' num2str(gamma)], 'Color', colors(idx));
    end

    xlabel('-log_{10}(\epsilon)');
    ylabel('Total SSE (Dual TLS)');
    title('Dual TLS Elbow Detection (Original Box-Counting)');
    legend('γ = 0.01', 'Elbow (γ = 0.01)', 'γ = 0.05', 'Elbow (γ = 0.05)');
    grid on; hold off;
end

function [m, b] = tls(x, y)
    x = x(:); y = y(:);
    X = [x - mean(x), y - mean(y)];
    [~, ~, V] = svd(X, 0);
    a = V(:, end);
    m = -a(1) / a(2);
    b = mean(y) - m * mean(x);
end
