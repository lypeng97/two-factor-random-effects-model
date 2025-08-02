function result = NR_optimization(epsilon, eta_ini, varpi_ini, max_iter, tol)
    % NR_OPTIMIZATION performs joint optimization of eta and varpi (interaction terms)
    %
    % Inputs:
    %   epsilon   - residual matrix (subjects × edges)
    %   eta_ini   - initial 1 × edges vector
    %   varpi_ini - initial subjects × 1 vector
    %   max_iter  - max number of iterations
    %   tol       - convergence tolerance
    %
    % Output:
    %   result    - struct containing final_eta, final_varpi, and objective_values

    % Initialize variables
    eta = eta_ini;                   % Initial eta values
    varpi = varpi_ini;               % Initial varpi values
    objective_values = zeros(max_iter, 1);  % Store objective values at each iteration
    objective_values(1) = sum((epsilon - varpi_ini * eta_ini) .^ 2, 'all');

    for iter = 2:max_iter
        % Step 1: Compute epsilon_new
        epsilon_new = varpi * eta;

        % Step 2: Compute gradients
        grad_eta = -2 * sum((epsilon - epsilon_new) .* varpi, 1);  % Gradient for eta
        grad_varpi = -2 * sum((epsilon - epsilon_new) .* eta, 2);  % Gradient for varpi
        grad = [grad_eta grad_varpi'];

        % Step 3: Compute Hessians
        s = size(varpi, 1);  % number of rows (subjects)
        I = size(eta, 2);    % number of columns (edges)

        hessian_eta = sum(2 * (varpi .^ 2)) * eye(I);
        hessian_varpi = sum(2 * (eta .^ 2)) * eye(s);
        hessian_cross = -2 * epsilon + 4 * epsilon_new;

        % Step 4: Combine Hessians into a full matrix
        H_upper = [hessian_eta, hessian_cross'];
        H_lower = [hessian_cross, hessian_varpi];
        full_H = [H_upper; H_lower];

        % Step 5: Invert Hessian if it's not too ill-conditioned
        if rcond(full_H) > 1e-16
            H_inv = inv(full_H);
        else
            disp('Stopping: Hessian is nearly singular.');
            break;
        end

        % Step 6: Update eta and varpi
        update_eta = grad * H_inv(1:I, :)';
        update_varpi = H_inv(I+1:end, :) * grad';

        eta = eta - update_eta;
        varpi = varpi - update_varpi;

        % Step 7: Update objective function value
        objective_values(iter) = sum((epsilon - varpi * eta) .^ 2, 'all');

        % Print progress
        fprintf('Iteration %d, Objective Value: %.6f\n', iter, objective_values(iter));

        % Step 8: Check convergence
        if norm(update_eta) < tol && norm(update_varpi) < tol
            fprintf('Convergence reached after %d iterations.\n', iter);
            break;
        end
    end

    % Return results
    result.final_eta = eta;
    result.final_varpi = varpi;
    result.objective_values = objective_values(1:iter);
end
