# Cart-Pendulum Optimal Control

Swing an inverted pendulum up from hanging and hold it upright. MATLAB + Simulink.

State: `x = [pos; vel; theta; thetadot]`, with `theta = 0` upright. Input: horizontal force on the cart.

## Stages

Built in three stages, each adding one thing:

1. **LQ**: infinite-horizon LQR, full state, holds near upright.
2. **LQ + EKF**: an extended Kalman filter estimates the state from a noisy pose (`pos`, `theta`), and the LQ gain runs on the estimate.
3. **Swing-up + LQ + EKF**: finite-horizon open-loop swing-up from hanging, then a supervisor hands over to LQ + EKF near the top.

## Files

| File | Role |
|------|------|
| `cartpole_dynamics.m` | Nonlinear dynamics `xdot = f(x,u)`. |
| `cartpole_setup.m` | Params, LQ design, swing-up lookup table. |
| `lqr_riccati.m` | Discrete LQR gain by backward Riccati. |
| `fh_controller.m` | Swing-up via multiple-shooting NLP (`fmincon`, SQP). |
| `swingup_lqr_fsm.m` | Supervisor: mode 0 swing-up, mode 1 LQR. |
| `compare_plan_vs_real.m` | Plan vs actual, tracking error. |
| Simulink models | One per stage. Hold the plant, EKF, and `CartPend` 3D view. |

## Params

```
M=0.5  m=0.2  l=0.3  b=0.1  g=9.81   [SI]
Ts=0.01        control sample time
Ts_swing=0.02  swing-up planning step
```

## Run

1. Run `cartpole_setup.m` (uncomment the finite-horizon block to build the swing-up table).
2. Open and run the Simulink model for the stage you want. Solver: `ode15s`, variable step.
3. Run `compare_plan_vs_real.m` after a swing-up run.

## Notes

- LQ is designed about `theta = 0`, so alone it only holds near the top. The swing-up handles the approach.
- Swing-up uses multiple shooting: states and inputs are both decision variables, dynamics enforced as equality constraints. Single shooting got stuck in local minima near upright; this fixed it. Still open loop, returns only the input.
- Supervisor catches when `|theta| < 0.15` and `|omega| < 1.5`, releases past `|theta| > 0.8`. The gap is hysteresis, so it does not chatter.
- Noise is added on the pose only. On for the EKF study, off for the swing-up run.

## Needs

MATLAB, Simulink, Optimization Toolbox, Control System Toolbox, Simscape Multibody.
