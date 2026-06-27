function [u_opt, x_pred] = fh_controller(x0, M, m, l, b, g, Ts_swing)
% FH_CONTROLLER  Finite-horizon swing-up via multiple-shooting NLP.
%   [u_opt, x_pred] = fh_controller(x0, M, m, l, b, g)
%
%   Solves an open-loop optimal control problem with fmincon (SQP),
%   optimizing inputs u_0..u_{N-1} and states x_1..x_N jointly.
%   RK4 dynamics defects enforced as equality constraints.
%
%   In  : x0          initial state [pos; vel; theta; omega]
%         M,m,l,b,g   plant parameters
%   Out : u_opt  N x 1     optimal input sequence
%         x_pred 4 x (N+1) optimal state trajectory
    x0 = x0(:);
    N     = 100;       % Horizon steps
    Ts    = Ts_swing; % Sampling time [s]
    sub   = 10;       % RK4 substeps per interval
    track = 0.8;      % cart position limit [m]
    V   = 2*diag([20 5 10 5]);          % stage state weight
    P   = 0.1;                          % input weight
    V_N = diag([5000, 500, 1000, 500]); % terminal weight
    nU = N;
    nX = 4*N;
% Initial guess: linear theta sweep to upright, rest zero (infeasible;
% equality constraints make it feasible).
    Xg      = zeros(4, N);
    Xg(3,:) = linspace(x0(3), 0, N);
    Ug      = zeros(N, 1);
    z0      = [Ug ; Xg(:)];
% Input saturation +/-20 N, states free
    lb = [-20*ones(nU,1) ; -inf(nX,1)];
    ub = [ 20*ones(nU,1) ;  inf(nX,1)];
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
'MaxIterations', 2000, 'MaxFunctionEvaluations', 3e5, ...
'StepTolerance', 1e-9, 'OptimalityTolerance', 1e-6, ...
'ConstraintTolerance', 1e-8);
    [z_opt, cost] = fmincon(@(z) ms_cost(z, x0, N, V, P, V_N), z0, ...
        [], [], [], [], lb, ub, ...
        @(z) ms_con(z, x0, N, Ts, sub, track, M, m, l, b, g), options); %#ok<ASGLU>
    u_opt  = z_opt(1:N);
    Xs     = reshape(z_opt(N+1:end), 4, N);
    x_pred = [x0, Xs];
    figure;
    t = 0:Ts:N*Ts;
    subplot(5,1,1); plot(t, x_pred(1,:), 'b', 'LineWidth', 1.5); ylabel('Cart Pos (m)'); grid on;
    subplot(5,1,2); plot(t, x_pred(2,:), 'r', 'LineWidth', 1.5); ylabel('Vel'); grid on;
    subplot(5,1,3); plot(t, x_pred(3,:), 'b', 'LineWidth', 1.5); ylabel('Theta (rad)'); grid on;
    subplot(5,1,4); plot(t, x_pred(4,:), 'r', 'LineWidth', 1.5); ylabel('Omega'); grid on;
    subplot(5,1,5); plot(t(1:end-1), u_opt, 'g', 'LineWidth', 1.5); ylabel('Input'); xlabel('Time (s)'); grid on;
    sgtitle('Open-Loop Optimal Trajectory (multiple shooting)');
end
function J = ms_cost(z, x0, N, V, P, V_N)
% Quadratic cost: sum(x'Vx + u'Pu) + x_N' V_N x_N
    U  = z(1:N);
    Xs = reshape(z(N+1:end), 4, N);
    J  = 0;
    x  = x0;
for i = 1:N
        J = J + (x' * V * x) + (U(i)' * P * U(i));
        x = Xs(:,i);
end
    J = J + (x' * V_N * x);
end
function [c, ceq] = ms_con(z, x0, N, Ts, sub, track, M, m, l, b, g)
% ceq: RK4 dynamics defects (4N).  c: cart track limits |pos| <= track (2N).
    U  = z(1:N);
    Xs = reshape(z(N+1:end), 4, N);
    ceq = zeros(4*N, 1);
    xp  = x0;
for i = 1:N
        xn = step_rk4(xp, U(i), Ts, sub, M, m, l, b, g);
        ceq((i-1)*4 + (1:4)) = Xs(:,i) - xn;
        xp = Xs(:,i);
end
    c = zeros(2*N, 1);
for i = 1:N
        c(2*i-1) =  Xs(1,i) - track;
        c(2*i)   = -Xs(1,i) - track;
end
end
function xn = step_rk4(x, u, Ts, sub, M, m, l, b, g)
% Fixed-step RK4 over one interval with `sub` substeps.
    h  = Ts/sub;
    xn = x;
for j = 1:sub
        k1 = cartpole_dynamics(xn,u,M,m,l,b,g);
        k2 = cartpole_dynamics(xn+h/2*k1,u,M,m,l,b,g);
        k3 = cartpole_dynamics(xn+h/2*k2,u,M,m,l,b,g);
        k4 = cartpole_dynamics(xn+h*k3,u,M,m,l,b,g);
        xn = xn + h/6*(k1 + 2*k2 + 2*k3 + k4);
end
end