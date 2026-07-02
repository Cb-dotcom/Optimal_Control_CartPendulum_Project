%% cartpole_setup.m
clear; clc;

%% Parameter Definition
M = 0.5;     % cart mass            [kg]
m = 0.2;     % pendulum mass        [kg]
l = 0.3;     % link length          [m]
b = 0.1;     % cart friction        [N/(m/s)]
g = 9.81;   % gravity              [m/s^2]
Ts = 0.01;     % sample time          [s]
x0 = [0; 0; pi; 0]; % initial state

%% Plant & Discretiziation
% continuous  xdot = f(x,u)
fc = @(x,u) cartpole_dynamics(x,u,M,m,l,b,g);

%% LQ
% Static continious linearization (Ac, Bc) about theta = 0
Ac = [ 0   1                 0                 0;
       0  -b/M              -m*g/M             0;
       0   0                 0                 1;
       0   b/(M*l)           (M+m)*g/(M*l)     0 ];
Bc = [ 0;
       1/M;
       0;
      -1/(M*l) ];
% Output matrix -> Measure the pose and  the EKF estimates the velocities.
G = [1 0 0 0;
     0 0 1 0];
% Discretize (ZOH)
% x_{k+1} = F x_k + W u_k
sys_d = c2d(ss(Ac,Bc,eye(4),0), Ts, 'zoh');
F = sys_d.A;    % discrete state matrix
W = sys_d.B;    % discrete input matrix
% Inf Horizon LQ regulator: constant gain Lbar (T_k -> Tbar)
% Designed offline- > in Simulink it is just one Gain block:  u = -L x.
V = diag([10 5 100 5]);       % state weight
P = 0.05;                     % input weight
L = lqr_riccati(F, W, V, P, 1e-12, 10000);   % = Lbar,  u = -L x

% %% FINITE HORIZON SET UP
% Ts_swing = 0.02;
% [u_opt, x_pred] = fh_controller(x0, M, m, l, b, g, Ts_swing);
% t_vec = (0:length(u_opt)-1)' * Ts_swing;

%% iLQR CONTROLLER OPEN LOOP SET UP
Ts_swing = 0.02;
[u_opt, x_pred] = iLQR_controller(x0, M, m, l, b, g, Ts_swing);
t_vec = (0:length(u_opt)-1)' * Ts_swing;

% Save it to the workspace
assignin('base', 't_swing_lut', t_vec);
assignin('base', 'u_swing_lut', u_opt);

%% Save to use in simulink
save('cartpole_cfg.mat','M','m','l','b','g','G','Ts','Ts_swing','fc','Ac','Bc','F','W','G','L','x0');
disp('Setup complete -> cartpole_cfg.mat saved.');