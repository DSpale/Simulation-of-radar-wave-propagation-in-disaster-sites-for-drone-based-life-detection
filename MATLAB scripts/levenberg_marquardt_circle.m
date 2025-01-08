function [center_est, radius_est] = levenberg_marquardt_circle(x_data, y_data, max_iter, tolerance, lambda)
    % Maximum number of iterations
    if ~exist('max_iter','var')
        max_iter = 100;
    end
    % Convergance tolerance
    if ~exist('tolerance','var')
        tolerance = 1e-6;
    end
    % Initial damping factor
    if ~exist('lambda','var')
        lambda = 1e-3;
    end

    % Initial guess for center and radius
    params = [mean(x_data), mean(y_data), std(x_data)]; % [cx, cy, r]
    num_points = length(x_data);
    if num_points ~= length(y_data)
        error('Number of x data and y data does not match up!')
    end
    
    % Function to compute residuals and Jacobian
    residuals = @(params) sqrt((x_data - params(1)).^2 + (y_data - params(2)).^2) - params(3);
    jacobian = @(params) [ % Analytical Jacobian matrix
        (params(1) - x_data)' ./ sqrt((params(1) - x_data).^2 + (params(2) - y_data).^2)', ...
        (params(2) - y_data)' ./ sqrt((params(1) - x_data).^2 + (params(2) - y_data).^2)', ...
        -ones(num_points, 1)
    ];
    
    % Iterative optimization loop
    for iter = 1:max_iter
        % Compute current residuals and Jacobian
        res = residuals(params)';
        J = jacobian(params);
    
        % Compute the normal equation and Levenberg-Marquardt step
        H = J' * J + lambda * eye(3); % Hessian approximation with damping
        g = J' * res; % Gradient
        dp = -H \ g; % Parameter update step
    
        % Update parameters and check for convergence
        params_new = params + dp';
        res_new = residuals(params_new)';
    
        % Check for improvement
        if norm(res_new) < norm(res)
            lambda = lambda / 10; % Decrease damping factor
            params = params_new;
            if norm(dp) < tolerance
                break; % Convergence achieved
            end
        else
            lambda = lambda * 10; % Increase damping factor
        end
    end
    % Extract the fitted parameters
    center_est = params(1:2);
    radius_est = params(3);
end