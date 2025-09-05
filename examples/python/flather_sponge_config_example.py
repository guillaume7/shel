import numpy as np


def build_flather_east_with_sponge_config(
    ny: int,
    nx: int,
    *,
    eta_east_value: float = 0.1,
    eta_relax: float = 0.5,
    momentum_relax: float = 0.3,
    sponge_width: int = 4,
    sponge_alpha: float = 0.2,
):
    """
    Build a conservative config combining Flather OBC (east) with a cosine-tapered sponge.

    Notes
    - Flather on the east boundary for both momentum and water level.
    - Conservative blending: eta_relax ~ 0.5 and momentum_relax ~ 0.3.
    - Sponge uses a mild alpha and cosine taper over a few cells.
    """
    eta_ext = np.zeros((ny, nx))
    eta_ext[:, -1] = eta_east_value
    return {
        "boundary_conditions": {
            "west": "closed",
            "east": "flather",
            "south": "closed",
            "north": "closed",
        },
        # External elevation specified on east boundary only
        "boundary_eta_ext": {"east": eta_ext},
        # Apply after continuity (default); set to "pre" to impose before the step
        "eta_bc_stage": "post",
        # Conservative relaxation toward external eta and momentum targets
        "boundary_eta_relax": float(eta_relax),
        "boundary_momentum_relax": float(momentum_relax),
        # Cosine-tapered sponge layer to smooth transition near the OBC
        "sponge": {
            "enabled": True,
            "width": int(sponge_width),  # number of interior cells to blend
            "alpha": float(sponge_alpha),  # blend strength (0..1), conservative ~0.1-0.3
            "taper": "cosine",  # cosine taper by default
            "apply_to": "both",  # blend eta and momentum
        },
    }


if __name__ == "__main__":
    cfg = build_flather_east_with_sponge_config(8, 10)
    print("Example config keys:", list(cfg.keys()))
    print("boundary_conditions:", cfg["boundary_conditions"])
    print("eta_bc_stage:", cfg.get("eta_bc_stage"))
    print("boundary_eta_relax:", cfg.get("boundary_eta_relax"))
    print("boundary_momentum_relax:", cfg.get("boundary_momentum_relax"))
    print("sponge:", cfg.get("sponge"))
