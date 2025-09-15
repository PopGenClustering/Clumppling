import os
import numpy as np
import matplotlib.pyplot as plt
from clumppling.plot import parse_custom_cmap, plot_colorbar, plot_alignment, plot_memberships_list, plot_graph
from clumppling.utils import str_to_pattern

K_max = 5
K_range = [3,5]
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

# Load membership matrices
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

Q_list_list = []
mode_names_list = []
for K in K_range:
    Q_list_list.append([Q for i, Q in enumerate(Q_list) if mode_K[i] == K])
    mode_names_list.append([name for i, name in enumerate(mode_names) if mode_K[i] == K])

# Load alignment results
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
cost_acrossK = {}
with open(os.path.join(acrossK_dir, "alignment_acrossK.txt"), "r") as f:
    next(f)
    for line in f:
        parts = line.strip().split(",")
        if len(parts) < 3:
            continue
        mode_pair = parts[0]
        alignment = str_to_pattern(parts[2])
        alignment_acrossK[mode_pair] = alignment
        cost = float(parts[1])
        cost_acrossK[mode_pair] = cost
        
# Plot alignment pattern
fig = plot_alignment(mode_K, mode_names, cmap, alignment_acrossK, all_modes_alignment, marker_size=200)
fig.savefig(os.path.join(fig_dir,"alignment_pattern_{}.png".format(suffix)), bbox_inches='tight', dpi=150, transparent=False)
plt.close(fig)

# Plot aligned memberships (as a list, with individuals ordered by major cluster's memberships, based on major mode of smallest K)
Q_ref = Q_list_list[0][0] # use the major mode of smallest K as reference
fig = plot_memberships_list(Q_list, cmap, names=mode_names, ind_labels=ind_labels, 
                            order_refQ=Q_ref, order_cls_by_label=True)
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, "aligned_membership_list.png"), bbox_inches='tight', dpi=150, transparent=False)
plt.close(fig)

# Plot aligned memberships (as a graph, with individuals ordered by major cluster's memberships, based on major mode of largest K)
Q_ref = Q_list_list[-1][0] # use the major mode of largest K as reference
fig = plot_graph(K_range, Q_list_list, cmap, 
                 names_list=mode_names_list, labels_list=None,
                 cost_acrossK=cost_acrossK, ind_labels=ind_labels, 
                 fontsize=14, line_cmap=plt.get_cmap("Greys"),
                 order_refQ=Q_ref, order_cls_by_label=True)
fig.savefig(os.path.join(fig_dir, "aligned_membership_graph.png"), bbox_inches='tight', dpi=150, transparent=False)
plt.close(fig)

# Plot aligned memberships (as a graph, with individuals unordered)
fig = plot_graph(K_range, Q_list_list, cmap, 
                 names_list=mode_names_list, labels_list=None,
                 cost_acrossK=cost_acrossK, ind_labels=ind_labels, 
                 fontsize=14, line_cmap=plt.get_cmap("Greys"),
                 order_refQ=None, order_cls_by_label=False)
fig.savefig(os.path.join(fig_dir, "aligned_membership_graph_unordered.png"), bbox_inches='tight', dpi=150, transparent=False)
plt.close(fig)
