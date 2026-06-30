%% compare_plan_vs_real.m
%  Overlay planned vs actual and tracking error.
%  Error is computed only during swing-up (mode == 0), since after the
%  catch the LQR drives and the open-loop plan no longer applies.

%% To ws variable names
sim_pos      = out.x;
sim_vel      = out.xdot;
sim_theta    = out.theta;
sim_omega    = out.thetadot;
sim_mode     = out.mode;

%% Extract time and data (each timeseries keeps its own time)
t_pos   = sim_pos.Time;      y_pos   = sim_pos.Data(:);
t_vel   = sim_vel.Time;      y_vel   = sim_vel.Data(:);
t_theta = sim_theta.Time;    y_theta = sim_theta.Data(:);
t_omega = sim_omega.Time;    y_omega = sim_omega.Data(:);
t_mode  = sim_mode.Time;     y_mode  = sim_mode.Data(:);

%% Planned trajectory
N  = length(u_opt);
Ts = Ts_swing;
t_plan = (0:N) * Ts;
T_plan_end = N * Ts;
p_pos   = x_pred(1,:);
p_vel   = x_pred(2,:);
p_theta = x_pred(3,:);
p_omega = x_pred(4,:);

%% Find when LQR catches (mode goes to 1)
idx_catch = find(y_mode >= 0.5, 1, 'first');
if ~isempty(idx_catch)
    t_catch = t_mode(idx_catch);
    fprintf('LQR catch at t = %.4f s\n', t_catch);
else
    t_catch = inf;
    fprintf('LQR never caught (mode stayed 0).\n');
end

% Compare only while swing-up is active: t <= min(t_catch, T_plan_end)
T_cmp = min(t_catch, T_plan_end);

%% Interpolate plan onto each measurement's time grid, up to T_cmp
mask_p = t_pos   <= T_cmp;   tc_p = t_pos(mask_p);
mask_v = t_vel   <= T_cmp;   tc_v = t_vel(mask_v);
mask_t = t_theta <= T_cmp;   tc_t = t_theta(mask_t);
mask_o = t_omega <= T_cmp;   tc_o = t_omega(mask_o);

e_pos   = y_pos(mask_p)   - interp1(t_plan, p_pos,   tc_p, 'linear', 'extrap');
e_vel   = y_vel(mask_v)   - interp1(t_plan, p_vel,   tc_v, 'linear', 'extrap');
e_theta = y_theta(mask_t) - interp1(t_plan, p_theta, tc_t, 'linear', 'extrap');
e_omega = y_omega(mask_o) - interp1(t_plan, p_omega, tc_o, 'linear', 'extrap');

%% Figs
labels = {'Pos [m]', 'Vel [m/s]', '\theta [rad]', '\omega [rad/s]'};
t_plans = {t_plan, t_plan, t_plan, t_plan};
p_data  = {p_pos, p_vel, p_theta, p_omega};
t_sims  = {t_pos, t_vel, t_theta, t_omega};
y_data  = {y_pos, y_vel, y_theta, y_omega};
t_errs  = {tc_p, tc_v, tc_t, tc_o};
e_data  = {e_pos, e_vel, e_theta, e_omega};

figure('Name','Plan vs Actual + Error','NumberTitle','off','Position',[60 60 1200 700]);

for k = 1:4
    % Left column: overlay (full sim time)
    subplot(5, 2, 2*k-1);
    plot(t_plans{k}, p_data{k}, 'b--', 'LineWidth', 1.8); hold on;
    plot(t_sims{k},  y_data{k}, 'r',   'LineWidth', 1.2);
    if isfinite(t_catch)
        xline(t_catch, 'k--', 'LineWidth', 1, 'Alpha', 0.6);
    end
    ylabel(labels{k}); grid on;
    if k == 1, title('Planned (blue) vs Actual (red)'); end
    if k == 1, legend('Planned','Actual','Catch','Location','best'); end
    if k == 4, xlabel('Time [s]'); end

    % Right column: error (swing-up phase only)
    subplot(5, 2, 2*k);
    plot(t_errs{k}, e_data{k}, 'm', 'LineWidth', 1.3);
    ylabel(['\Delta ' labels{k}]); grid on;
    if k == 1, title(sprintf('Tracking Error (swing-up only, 0 to %.2f s)', T_cmp)); end
    if k == 4, xlabel('Time [s]'); end
end

% Bottom row: mode
subplot(5, 2, [9 10]);
stairs(t_mode, y_mode, 'g', 'LineWidth', 1.5);
if isfinite(t_catch)
    xline(t_catch, 'k--', 'LineWidth', 1, 'Alpha', 0.6);
end
ylabel('Mode'); xlabel('Time [s]');
yticks([0 1]); yticklabels({'Swing-up','LQR'}); grid on;
ylim([-0.2 1.3]);

%% Print summary
fprintf('\n=== Tracking Error (swing-up phase: 0 to %.3f s) ===\n', T_cmp);
fprintf('  pos   :  max|e| = %.4f m      rms = %.4f m\n',     max(abs(e_pos)),   rms(e_pos));
fprintf('  vel   :  max|e| = %.4f m/s    rms = %.4f m/s\n',   max(abs(e_vel)),   rms(e_vel));
fprintf('  theta :  max|e| = %.4f rad    rms = %.4f rad\n',   max(abs(e_theta)), rms(e_theta));
fprintf('  omega :  max|e| = %.4f rad/s  rms = %.4f rad/s\n', max(abs(e_omega)), rms(e_omega));
fprintf('  LQR catch time : %.3f s  (plan horizon: %.2f s)\n', t_catch, T_plan_end);
fprintf('===================================================\n');