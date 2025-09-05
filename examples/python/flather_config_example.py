import numpy as np


def build_flather_east_config(ny: int, nx: int, eta_east_value: float = 0.1):
    eta_ext = np.zeros((ny, nx))
    eta_ext[:, -1] = eta_east_value
    return {
        "boundary_conditions": {"west": "closed", "east": "flather", "south": "closed", "north": "closed"},
        "boundary_eta_ext": {"east": eta_ext},
        # Apply after continuity (default); set to "pre" to impose before the step
        "eta_bc_stage": "post",
        # 1.0 = pure Dirichlet at eta boundary; <1.0 = relaxation toward external value
        "boundary_eta_relax": 1.0,
    }


if __name__ == "__main__":
    cfg = build_flather_east_config(8, 10)
    print("Example config keys:", list(cfg.keys()))
    for k, v in cfg["boundary_conditions"].items():
        print(f"  {k}: {v}")
