import os
import numpy as np
import matplotlib.pyplot as plt
from clumppling.plot import parse_custom_cmap, plot_colorbar, plot_alignment, plot_memberships_list
from clumppling.utils import str_to_pattern

K_max = 5
cmap_file = "examples/custom_colors.txt"
acrossK_dir = "examples/submodules/K3K5_acrossK_output"
mode_dir = os.path.join(acrossK_dir, "modes_aligned")
ind_label_file = "examples/capeverde_ind_labels.txt"
fig_dir = "examples/submodules/figures"
os.makedirs(fig_dir, exist_ok=True)

# Read the individual labels
ind_labels = []
with open(ind_label_file, "r") as f:
    for line in f:
        ind_labels.append(line.strip())
print(f"Loaded {len(ind_labels)} individual labels from {ind_label_file}")

# Read the custom color map from the file
with open(cmap_file, "r") as f:
    colors = f.read().strip().splitlines()
cmap = parse_custom_cmap(colors, K=K_max)
# Plot the colorbar
plot_colorbar(cmap,K_max,fig_dir)

# plot alignment pattern
if not os.path.exists(mode_dir):
    raise FileNotFoundError(f"Mode directory {mode_dir} does not exist.")
Q_list = []
mode_names = []
for f in os.listdir(mode_dir):
    if f.endswith(".Q"):
        suffix = f.strip(".Q")
        Q_list.append(np.loadtxt(os.path.join(mode_dir, f)))
        mode_names.append(suffix)
mode_K = [Q.shape[1] for Q in Q_list]
assert mode_K == sorted(mode_K), "Modes should be sorted by K."

all_modes_alignment = {}
with open(os.path.join(mode_dir, "all_modes_alignment.txt"), "r") as f:
    for line in f:
        parts = line.strip().split(":")
        if len(parts) < 2:
            continue
        mode_name = parts[0]
        alignment = str_to_pattern(parts[1])
        all_modes_alignment[mode_name] = alignment

alignment_acrossK = {}
with open(os.path.join(acrossK_dir, "alignment_acrossK.txt"), "r") as f:
    next(f)
    for line in f:
        parts = line.strip().split(",")
        if len(parts) < 3:
            continue
        mode_pair = parts[0]
        alignment = str_to_pattern(parts[2])
        alignment_acrossK[mode_pair] = alignment

fig = plot_alignment(mode_K, mode_names, cmap, alignment_acrossK, all_modes_alignment, marker_size=200)
fig.savefig(os.path.join(fig_dir,"alignment_pattern_{}.png".format(suffix)), bbox_inches='tight', dpi=150, transparent=False)
plt.close(fig)

# plot aligned memberships (as a list)
Q_ref = Q_list[-1]
fig = plot_memberships_list(Q_list, cmap, names=mode_names, ind_labels=ind_labels, 
                            order_refQ=Q_ref, order_cls_by_label=True)
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, "aligned_membership.png"), bbox_inches='tight', dpi=150, transparent=False)
plt.close(fig)
