import numpy as np

from shel.model.boundary_conditions import get_bc


def test_sommerfeld_eta_formula_east_west():
    ny, nx = 4, 8
    H = np.ones((ny, nx)) * 16.0
    g = 9.81
    dt = 0.1
    dx = 2.0
    dy = 3.0

    e_cls = get_bc("radiative")[1]
    assert e_cls is not None
    bc = e_cls()

    # Construct a linear field so the one-sided difference is constant
    x = np.arange(nx)
    eta_old = np.tile(x, (ny, 1)).astype(float)
    eta_next = eta_old.copy()

    c = np.sqrt(g * H.mean())
    r = c * dt / dx

    bc.apply_side_eta(eta_next, "west", eta_old=eta_old, H=H, g=g, dt=dt, dx=dx, dy=dy)
    assert np.allclose(
        eta_next[:, 0], eta_old[:, 0] - r * (eta_old[:, 1] - eta_old[:, 0])
    )

    eta_next2 = eta_old.copy()
    bc.apply_side_eta(eta_next2, "east", eta_old=eta_old, H=H, g=g, dt=dt, dx=dx, dy=dy)
    assert np.allclose(
        eta_next2[:, -1], eta_old[:, -1] - r * (eta_old[:, -1] - eta_old[:, -2])
    )


def test_sommerfeld_eta_formula_south_north():
    ny, nx = 6, 5
    H = np.ones((ny, nx)) * 9.0
    g = 9.81
    dt = 0.05
    dx = 1.5
    dy = 2.0

    e_cls = get_bc("radiative")[1]
    assert e_cls is not None
    bc = e_cls()

    y = np.arange(ny)
    eta_old = np.tile(y.reshape(-1, 1), (1, nx)).astype(float)
    eta_next = eta_old.copy()

    c = np.sqrt(g * H.mean())
    r = c * dt / dy

    bc.apply_side_eta(eta_next, "south", eta_old=eta_old, H=H, g=g, dt=dt, dx=dx, dy=dy)
    assert np.allclose(
        eta_next[0, :], eta_old[0, :] - r * (eta_old[1, :] - eta_old[0, :])
    )

    eta_next2 = eta_old.copy()
    bc.apply_side_eta(
        eta_next2, "north", eta_old=eta_old, H=H, g=g, dt=dt, dx=dx, dy=dy
    )
    assert np.allclose(
        eta_next2[-1, :], eta_old[-1, :] - r * (eta_old[-1, :] - eta_old[-2, :])
    )


def test_sommerfeld_pulse_propagation_right_going():
    # 1D-like right-going pulse should advance eastward roughly c*dt/dx per step under linearized SWE
    ny, nx = 8, 120
    dx = dy = 1.0
    H = np.ones((ny, nx)) * 10.0
    g = 9.81
    c = np.sqrt(g * H.mean())
    dt = 0.05  # r ~ 0.495
    r = c * dt / dx

    # Initial Gaussian eta centered away from boundaries
    x = np.arange(nx)
    X = np.tile(x, (ny, 1))
    x0 = int(nx * 0.33)
    eta = 0.02 * np.exp(-((X - x0) ** 2) / (2.0 * 5.0**2))

    # Initialize a right-going mode: U ~ (c/H) * eta averaged to faces; V = 0
    U = np.zeros((ny, nx + 1))
    V = np.zeros((ny + 1, nx))
    eta_face = 0.5 * (eta[:, :-1] + eta[:, 1:])
    U[:, 1:-1] = (c / H.mean()) * eta_face

    from shel.model.solvers.common.stepper import explicit_step_with_config

    cfg = {
        "boundary_conditions": {
            "west": "radiative",
            "east": "radiative",
            "south": "closed",
            "north": "closed",
        }
    }

    # Track centroid of positive eta along x
    def centroid_x(field):
        prof = field.mean(axis=0)
        prof = np.where(prof > 0, prof, 0.0)
        s = prof.sum()
        if s <= 0:
            return float(np.argmax(prof))
        return float((np.arange(field.shape[1]) * prof).sum() / s)

    c0 = centroid_x(eta)
    steps = 40
    for _ in range(steps):
        eta, U, V, _ = explicit_step_with_config(
            eta,
            H,
            U,
            V,
            dt=dt,
            dx=dx,
            dy=dy,
            g=g,
            r=0.0,
            nu=0.0,
            enable_advection=False,
            f=None,
            enable_coriolis=False,
            config=cfg,
        )
    c1 = centroid_x(eta)

    expected_shift = steps * r
    measured_shift = c1 - c0
    # Loose tolerance: numerical dispersion + simple initializer
    assert measured_shift > 0
    assert abs(measured_shift - expected_shift) < 3.0  # within ~3 cells
