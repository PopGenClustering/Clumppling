import os
import numpy as np
import matplotlib.pyplot as plt
from clumppling.plot import parse_custom_cmap,plot_membership, plot_memberships_list

runs_dir = "examples/capeverde_output/input"
cmap_file = "clumppling/files/default_colors.txt"
fig_dir = "examples/unaligned_runs/figures"
os.makedirs(fig_dir, exist_ok=True)

# load all Q matrices
Q_list = []
run_names = []
for f in os.listdir(runs_dir):
    if f.endswith(".Q"):
        suffix = f.strip(".Q")
        Q_list.append(np.loadtxt(os.path.join(runs_dir, f)))
        run_names.append(suffix)

K_max = max(Q.shape[1] for Q in Q_list)
# Read the custom color map from the file
with open(cmap_file, "r") as f:
    colors = f.read().strip().splitlines()
cmap = parse_custom_cmap(colors, K=K_max)

n_cols = 4
n_rows = int(np.ceil(len(Q_list) / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, n_rows))
axes = axes.flatten()

for i in range(len(Q_list)):
    Q = Q_list[i]
    ax = axes[i] if len(Q_list)>1 else axes
    plot_membership(Q, cmap, ax=ax, ylab="")
    
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, "unaligned_runs_memberships.png"), dpi=300)
plt.close(fig)