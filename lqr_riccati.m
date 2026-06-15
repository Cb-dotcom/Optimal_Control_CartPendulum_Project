function [L, T] = lqr_riccati(F, W, V, P, tol, maxit)
%   LQR_RICCATI  Infinite-horizon discrete LQR
%
%
%
%     System : x_{k+1} = F x_k + W u_k
%     Cost J : sum_k ( x_k' V x_k + u_k' P u_k )
%     Gain   : u_k = -L x_k      (returned L = Lbar)
%     T      : converged Riccati matrix (= Tbar)
%
%   Recursion until we converge:
%     L_{k}  = (P + W'T_{k+1}W)^{-1} W'T_{k+1}F
%     T_{k}  =  V + F'T_{k+1}F - F'T_{k+1}W L_{k}

    % start from T_N = V 
    T = V;                                  
    for k = 1:maxit
        L    = (P + W'*T*W) \ (W'*T*F);     % gain using current T (= T_{k+1})
        Tnew = V + F'*T*F - F'*T*W*L;       % one step back -> T_k
        if norm(Tnew - T, 'fro') < tol
            T = Tnew;
            L = (P + W'*T*W) \ (W'*T*F);    % final gain from the converged T
            return;
        end
        T = Tnew;
    end
    warning('lqr_riccati:noConverge', ...
            'Riccati did not converge in %d iterations.', maxit);
    L = (P + W'*T*W) \ (W'*T*F);
end