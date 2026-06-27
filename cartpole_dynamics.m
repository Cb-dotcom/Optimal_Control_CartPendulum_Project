function xdot = cartpole_dynamics(x, u, M, m, l, b, g)
%   Continuous-time nonlinear cart-pole
%
%   xdot = cartpole_dynamics(x, u, M, m, l, b, g)
%
%
%   State : x = [ pos ; vel ; theta ; omega ]
%             pos: cart position        [m]
%             vel: cart velocity        [m/s]
%             theta: pendulum angle     [rad]   (theta = 0 is upright)
%             omega: pendulum angular vel      [rad/s]
%   Input : u = F:   horizontal force on cart  [N]
%   Params: M, m, l, b, g
%
    vel = x(2);
    th  = x(3);
    om  = x(4);
    F = u(1);
    s = sin(th);  c = cos(th);
    den    = M + m*s^2;
    common = F - b*vel + m*l*s*om^2;
% d(vel)/dt
    xacc  = ( common - m*g*c*s ) / den;
% d(omega)/dt
    thacc = ( (M+m)*g*s - c*common ) / (l*den);
    xdot = [ vel ; xacc ; om ; thacc ];
end