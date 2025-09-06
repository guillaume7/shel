import numpy as np

from shel.model.solvers.common.ministep import explicit_step
from shel.model.solvers.time import leapfrog_stepper


def test_leapfrog_inertial_oscillation_quarter_period():
    ny, nx = 40, 40
    dx = dy = 1.0
    dt = 0.05
    g = 9.81

    # Constant depth and f-plane
    H = np.ones((ny, nx)) * 10.0
    f = np.ones((ny, nx)) * 0.1  # rad/s

    # Initial state: eta=0, uniform interior U=1, V=0
    eta0 = np.zeros((ny, nx))
    U0 = np.zeros((ny, nx + 1))
    V0 = np.zeros((ny + 1, nx))
    U0[:, 1:-1] = 1.0

    # Startup: one explicit Euler step to obtain state at n from n-1
    eta1, U1, V1 = explicit_step(
        eta0,
        H,
        U0,
        V0,
        dt=dt,
        dx=dx,
        dy=dy,
        g=g,
        r=0.0,
        nu=0.0,
        enable_advection=False,
        f=f,
        enable_coriolis=True,
        bc_type="closed",
    )

    # Target quarter inertial period: T=2*pi/f; quarter at pi/2
    omega = f.mean()  # constant
    target_time = 0.5 * np.pi / omega
    nsteps = int(round(target_time / dt))

    eta_nm1, U_nm1, V_nm1 = eta0, U0, V0
    eta_n, U_n, V_n = eta1, U1, V1

    for _ in range(nsteps):
        eta_np1, U_np1, V_np1, eta_n_f, U_n_f, V_n_f = leapfrog_stepper(
            eta_nm1,
            eta_n,
            H,
            U_nm1,
            U_n,
            V_nm1,
            V_n,
            dt=dt,
            dx=dx,
            dy=dy,
            g=g,
            r=0.0,
            nu_visc=0.0,
            enable_advection=False,
            f=f,
            enable_coriolis=True,
            bc_type="closed",
            asselin_nu=0.02,
        )
        # advance time levels
        eta_nm1, U_nm1, V_nm1 = eta_n, U_n, V_n
        eta_n, U_n, V_n = eta_np1, U_np1, V_np1

    # Measure interior means (avoid boundary effects)
    u_mean = U_n[:, 1:-1].mean()
    v_mean = V_n[1:-1, :].mean()

    # Smoke-level checks: boundedness and rotational response (v becomes negative),
    # and some decay in u from the initial value due to rotation + filtering.
    assert np.isfinite(u_mean) and np.isfinite(v_mean)
    assert abs(u_mean) < 1.1 and abs(v_mean) < 1.1
    assert v_mean < -0.1
    assert u_mean < 0.9
