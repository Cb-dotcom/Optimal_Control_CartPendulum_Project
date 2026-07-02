function [u_opt, x_pred] = iLQR_controller(x0, M, m, l, b, g, Ts_swing)
% =========================================================================
% iLQR_controller - Iterative Linear Quadratic Regulator (Open Loop)
%
% INPUTS:
%   x0       - initial state [pos; vel; theta; omega]
%   M,m,l,b,g - cart-pole physical parameters
%   Ts_swing  - sample time [s]
%
% OUTPUTS:
%   u_opt  - optimal input sequence (N x 1)
%   x_pred - predicted state trajectory (4 x N+1)
%
% LOGIC:
%   iLQR solves a nonlinear optimal control problem via Dynamic Programming.
%   It iterates two passes until convergence:
%
%   1. BACKWARD PASS (Riccati):
%      Linearize dynamics around current trajectory at each step k.
%      Propagate cost-to-go backwards to compute:
%        - d_k : feedforward gain (drives system toward goal)
%        - K_k : feedback gain   (corrects trajectory deviations)
%
%   2. FORWARD PASS (Rollout):
%      Apply updated control law:
%        u_k = u_old_k + alpha*d_k + K_k*(x_k - x_old_k)
%      Simulate forward with true nonlinear dynamics (RK4).
%      Line search on alpha ensures cost decreases.
%
%   Unlike LQ (u = -Lx around equilibrium), iLQR linearizes around a
%   reference trajectory, so the gain corrects deviations from it.
%   As it converges: d_k -> 0 and x_k -> x_old_k.
% =========================================================================

    x0 = x0(:);
    N   = 100;
    Ts  = Ts_swing;
    sub = 10;

    % Weights
    V   = 2*diag([20 5 10 5]);
    P   = 0.1;
    V_N = diag([5000, 500, 1000, 500]);

    % Initial trajectory
    X_traj      = zeros(4, N+1);
    X_traj(:,1) = x0;
    U_traj      = zeros(N, 1);
    for i = 1:N
        X_traj(:,i+1) = step_rk4(X_traj(:,i), U_traj(i), Ts, sub, M, m, l, b, g);
    end

    J_old = compute_cost(X_traj, U_traj, N, V, P, V_N);

    % iLQR main loop
    for iter = 1:200

        % --- BACKWARD PASS ---
        T_k = V_N;
        t_k = V_N * X_traj(:, N+1);  % gradient of terminal cost

        K_gain = zeros(N, 4);
        d_gain = zeros(N, 1);

        for k = N:-1:1
            x_k    = X_traj(:, k);
            u_k    = U_traj(k);

            [F_k, W_k] = get_jacobians(x_k, u_k, Ts, sub, M, m, l, b, g);

            % Helper Matrices to compute the optimization feedback gain K and d
            Q_x = V * x_k + F_k' * t_k;
            Q_u = P * u_k + W_k' * t_k;
            Q_uu = P   + W_k' * T_k * W_k;
            Q_ux = W_k' * T_k * F_k;
            Q_xx = V   + F_k' * T_k * F_k;

            % Gains (matches: L = (P + W'TW)^{-1} W'TF)
            d_gain(k)   = -Q_uu \ Q_u;
            K_gain(k,:) = -Q_uu \ Q_ux;

            % Riccati update (matches: T_k = V + F'T F - F'TW L)
            T_k = Q_xx - K_gain(k,:)' * Q_uu * K_gain(k,:);
            t_k = Q_x  - K_gain(k,:)' * Q_uu * d_gain(k);
        end

        % --- FORWARD PASS ---
        alpha  = 1.0;
        X_new  = zeros(4, N+1);
        U_new  = zeros(N, 1);
        X_new(:,1) = x0;

        for k = 1:N
            %compute the optimal control action using the precomputed gains
            U_new(k)      = U_traj(k) + alpha * d_gain(k) + K_gain(k,:) * (X_new(:,k) - X_traj(:,k));
            %propagation of the state
            X_new(:,k+1)  = step_rk4(X_new(:,k), U_new(k), Ts, sub, M, m, l, b, g);
        end

        J_new = compute_cost(X_new, U_new, N, V, P, V_N);

        % Line search to keep optimize until convergence
        while J_new >= J_old && alpha > 1e-4
            alpha = alpha * 0.5;
            for k = 1:N
                U_new(k)     = U_traj(k) + alpha * d_gain(k) + K_gain(k,:) * (X_new(:,k) - X_traj(:,k));
                X_new(:,k+1) = step_rk4(X_new(:,k), U_new(k), Ts, sub, M, m, l, b, g);
            end
            J_new = compute_cost(X_new, U_new, N, V, P, V_N);
        end

        if J_new >= J_old
            fprintf('Converged at iteration %d\n', iter);
            break;
        end

        X_traj = X_new;
        U_traj = U_new;
        J_old  = J_new;
    end

    u_opt  = U_traj;
    x_pred = X_traj;

    % Plot
    figure;
    t = 0:Ts_swing:N*Ts_swing;
    subplot(5,1,1); plot(t, x_pred(1,:), 'b', 'LineWidth', 1.5); ylabel('Cart Pos (m)'); grid on;
    subplot(5,1,2); plot(t, x_pred(2,:), 'r', 'LineWidth', 1.5); ylabel('Vel');           grid on;
    subplot(5,1,3); plot(t, x_pred(3,:), 'b', 'LineWidth', 1.5); ylabel('Theta (rad)');   grid on;
    subplot(5,1,4); plot(t, x_pred(4,:), 'r', 'LineWidth', 1.5); ylabel('Omega');         grid on;
    subplot(5,1,5); plot(t(1:end-1), u_opt, 'g', 'LineWidth', 1.5); ylabel('Input'); xlabel('Time (s)'); grid on;
    sgtitle('Open-Loop Optimal Trajectory (iLQR)');
end

% --- HELPER FUNCTIONS ---

function J = compute_cost(X, U, N, V, P, V_N)
    J = 0;
    for i = 1:N
        J = J + X(:,i)'*V*X(:,i) + U(i)'*P*U(i);
    end
    J = J + X(:,N+1)'*V_N*X(:,N+1);
end
% Compute the linearization of the non linear dynamics
function [A, B] = get_jacobians(x, u, Ts, sub, M, m, l, b, g)
    eps = 1e-5;
    nx  = length(x);
    f0  = step_rk4(x, u, Ts, sub, M, m, l, b, g);

    A = zeros(nx, nx);
    for i = 1:nx
        x_eps    = x; x_eps(i) = x_eps(i) + eps;
        A(:,i)   = (step_rk4(x_eps, u, Ts, sub, M, m, l, b, g) - f0) / eps;
    end

    B = zeros(nx, 1);
    u_eps = u + eps;
    B     = (step_rk4(x, u_eps, Ts, sub, M, m, l, b, g) - f0) / eps;
end

function xn = step_rk4(x, u, Ts, sub, M, m, l, b, g)
    h  = Ts/sub;
    xn = x;
    for j = 1:sub
        k1 = cartpole_dynamics(xn,          u, M, m, l, b, g);
        k2 = cartpole_dynamics(xn+h/2*k1,   u, M, m, l, b, g);
        k3 = cartpole_dynamics(xn+h/2*k2,   u, M, m, l, b, g);
        k4 = cartpole_dynamics(xn+h*k3,     u, M, m, l, b, g);
        xn = xn + h/6*(k1 + 2*k2 + 2*k3 + k4);
    end
end

function dx = cartpole_dynamics(x, u, M, m, l, b, g)
    st = sin(x(3)); ct = cos(x(3));
    dxdt    = zeros(4,1);
    denom   = M + m*st^2;
    dxdt(1) = x(2);
    dxdt(2) = (m*l*x(4)^2*st - m*g*st*ct + u - b*x(2)) / denom;
    dxdt(3) = x(4);
    dxdt(4) = (-m*l*x(4)^2*st*ct + (M+m)*g*st - ct*(u - b*x(2))) / (l*denom);
    dx      = dxdt;
end