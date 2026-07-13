import cvxpy as cp
import numpy as np

from clumppling.core import align_ILP


def test_glpk_mi_is_available() -> None:
    """The solver required by Clumppling must be installed."""
    assert "GLPK_MI" in cp.installed_solvers()


def test_align_ilp_recovers_known_permutation() -> None:
    """Clumppling should recover a known swap of two clusters."""
    p = np.array(
        [
            [0.9, 0.1],
            [0.8, 0.2],
            [0.1, 0.9],
            [0.2, 0.8],
        ],
        dtype=float,
    )

    # Swap the two ancestry columns.
    q = p[:, [1, 0]]

    objective, q_to_p = align_ILP(p, q)

    assert objective is not None
    assert np.isclose(objective, 0.0, atol=1e-10)
    assert q_to_p == [1, 0]
