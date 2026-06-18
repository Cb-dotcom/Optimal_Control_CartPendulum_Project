function u_cmd = run_mpc_optimization(x_hat)
    % --- MPC Parameters ---
    N = 40; % Horizon length
    
    % Penalty Matrices
    V = diag([10, 1, 50, 1]); % State penalty (Heavy penalty on angle theta)
    P = 1;                    % Control effort penalty (u^2)
    V_N = diag([50, 10, 200, 10]); % Terminal cost h_N 
    
    % Persistent variable for Warm Starting
    persistent u_prev;
    if isempty(u_prev)
        u_prev = zeros(N, 1); % Initial guess is zero effort
    end
    
    % --- Warm Start Shifting ---
    u_guess = zeros(N, 1);
    u_guess(1:N-1) = u_prev(2:N); 
    u_guess(N) = u_prev(N); 
    
    % --- Solver Options ---
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'none');
    
    % --- Actuator Constraints ---
    u_min = -12 * ones(N, 1); % -12V limit
    u_max = 12 * ones(N, 1);  % +12V limit
    
    % --- Solve the Optimization ---
    % Because we are native MATLAB now, these function handles are perfectly legal!
    u_opt = fmincon(@(U) compute_cost(U, x_hat, N, V, P, V_N), ...
                    u_guess, [], [], [], [], u_min, u_max, ...
                    @(U) track_constraints(U, x_hat, N), options);
    
    % Save sequence for next step's warm start
    if ~isempty(u_opt)
        u_prev = u_opt;
        u_cmd = u_opt(1); % Receding Horizon: Apply only the first input
    else
        % Fallback if fmincon fails to find a solution
        u_cmd = 0; 
    end
end

% =========================================================================
% Helper Functions (Stored at the bottom of the same .m file)
% =========================================================================

function J = compute_cost(U, x_init, N, V, P, V_N)
    J = 0;
    x_current = x_init;
    Ts = 0.05; % Prediction step size
    
    for i = 1:N
        u_i = U(i);
        J = J + (x_current' * V * x_current) + (u_i' * P * u_i);
        
        x_dot = plant_dynamics(x_current, u_i);
        x_current = x_current + x_dot * Ts; % Euler integration
    end
    J = J + (x_current' * V_N * x_current);
end

function [c, ceq] = track_constraints(U, x_init, N)
    ceq = []; 
    c = zeros(2*N, 1); 
    
    x_current = x_init;
    Ts = 0.05;
    max_track_length = 0.8; 
    
    for i = 1:N
        u_i = U(i);
        x_dot = plant_dynamics(x_current, u_i);
        x_current = x_current + x_dot * Ts;
        
        c(2*i - 1) = x_current(1) - max_track_length; 
        c(2*i) = -x_current(1) - max_track_length;    
    end
end

function x_dot = plant_dynamics(x, u)
    % Ensure these parameters match your physical plant exactly
    M = 0.5; m = 0.2; l = 0.3; b = 0.1; g = 9.81;
    vel = x(2); th = x(3); om = x(4);
    
    s = sin(th); c = cos(th);
    den = M + m*s^2;
    common = u - b*vel + m*l*s*om^2;
    
    xacc = (common - m*g*c*s) / den;
    thacc = ((M+m)*g*s - c*common) / (l*den);
    
    x_dot = [vel; xacc; om; thacc];
end