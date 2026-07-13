import numpy as np

from clumppling.core import cd_default


def test_louvain_is_reproducible_with_fixed_seed() -> None:
    """Louvain should return the same partition with the same seed."""

    # Two strongly connected groups joined by one weak edge.
    adjacency = np.array(
        [
            [0.00, 1.00, 1.00, 0.01, 0.00, 0.00],
            [1.00, 0.00, 1.00, 0.00, 0.00, 0.00],
            [1.00, 1.00, 0.00, 0.00, 0.00, 0.00],
            [0.01, 0.00, 0.00, 0.00, 1.00, 1.00],
            [0.00, 0.00, 0.00, 1.00, 0.00, 1.00],
            [0.00, 0.00, 0.00, 1.00, 1.00, 0.00],
        ],
        dtype=float,
    )

    np.random.seed(42)
    partition_1 = cd_default(
        adjacency,
        method="louvain",
        res=1.0,
    )

    np.random.seed(42)
    partition_2 = cd_default(
        adjacency,
        method="louvain",
        res=1.0,
    )

    assert partition_1 == partition_2
    assert len(set(partition_1)) == 2
