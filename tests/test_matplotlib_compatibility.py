import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from clumppling.plot import plot_graph


def test_plot_graph_default_greys_colormap():
    """Default connection-line colors work with supported Matplotlib versions."""

    q_k2 = np.array([
        [0.6, 0.4],
        [0.3, 0.7],
    ])

    q_k3_mode1 = np.array([
        [0.5, 0.3, 0.2],
        [0.2, 0.3, 0.5],
    ])

    q_k3_mode2 = np.array([
        [0.4, 0.4, 0.2],
        [0.1, 0.2, 0.7],
    ])

    fig = plot_graph(
        K_range=[2, 3],
        Q_list_list=[
            [q_k2],
            [q_k3_mode1, q_k3_mode2],
        ],
        cmap=[
            (0.8, 0.1, 0.1),
            (0.1, 0.8, 0.1),
            (0.1, 0.1, 0.8),
        ],
        names_list=[
            ["K2M1"],
            ["K3M1", "K3M2"],
        ],
        cost_acrossK={
            "K2M1-K3M1": 0.1,
            "K2M1-K3M2": 0.2,
        },
        alt_color=False,
        order_cls_by_label=False,
    )

    try:
        assert fig is not None
        assert len(fig.axes) == 3
    finally:
        plt.close(fig)
