function xdot = cartpole_dynamics(x, u, p)
%   Continuous-time nonlinear cart-pole
%
%   xdot = cartpole_dynamics(x, u, p)
%
%
%   State : x = [ pos ; vel ; theta ; omega ]   
%             pos: cart position        [m]
%             vel: cart velocity        [m/s]
%             theta: pendulum angle     [rad]   (theta = 0 is upright)
%             omega: pendulum angular vel      [rad/s]
%   Input : u = F:   horizontal force on cart  [N]
%   Params: p.M, p.m, p.l, p.b, p.g
%

    vel = x(2);
    th  = x(3);
    om  = x(4);

    M = p.M;  m = p.m;  l = p.l;  b = p.b;  g = p.g;
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