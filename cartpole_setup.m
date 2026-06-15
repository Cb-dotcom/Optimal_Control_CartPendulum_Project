%% cartpole_setup.m
% Methodology:
%   1) True Plant  = continuous nonlinear cart-pole (cartpole_dynamics.m)
%   2) Discritize the true plant using forward euler (Ts in the 10-40 ms range, avoid 5)
%   3) Linearize and use LQ, then add noise  and use EKF (check LQ+KF vs
%   LQ+EKF)
%   4) Solve finite-horizon / moving-horizon optimal control on the
%      discitized non linear -> NLQ and MPC + EKF when we add noise.

clear; clc;

%%
p.M = 0.5;     % cart mass            [kg]
p.m = 0.2;     % pendulum mass        [kg]
p.l = 0.3;     % link length          [m]
p.b = 0.1;     % cart friction        [N/(m/s)]
p.g = 9.81;   % gravity              [m/s^2]

Ts = 0.01;     % sample time          [s]

%% Plant & Discretiziation 
% continuous  xdot = f(x,u)
fc = @(x,u) cartpole_dynamics(x,u,p);        
% discrete (forward Euler)
fd = @(x,u) x + Ts*fc(x,u);                   


%% LQ 
% Static continious linearization (Ac, Bc) about theta = 0
Ac = [ 0   1                       0                          0;    
       0  -p.b/p.M                -p.m*p.g/p.M                0;     
       0   0                       0                          1;     
       0   p.b/(p.M*p.l)           (p.M+p.m)*p.g/(p.M*p.l)    0 ];   
 
Bc = [ 0;
       1/p.M;           
       0;
      -1/(p.M*p.l) ];    
 
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
V = diag([10 5 100 5]);       %state weight
P = 0.05;                     % input weight 
L = lqr_riccati(F, W, V, P, 1e-12, 10000);   % = Lbar,  u = -L x


%% Save to use in simulink
save('cartpole_cfg.mat','p','Ts','fc','fd','Ac','Bc','F','W','G','L');
disp('Setup complete -> cartpole_cfg.mat saved.');
