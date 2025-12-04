# !/bin/bash

# default run with default: louvain mode detection, representative Q matrices for modes, and alignment across-K using the best mode pair
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_output \
-f admixture \
--extension .indivq \
--ind_labels examples/capeverde_ind_labels.txt

# run with average Q matrices for modes and alignment across-K using the major mode pair
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_avg_output \
-f admixture \
--use_rep F \
--use_best_pair F \
--extension .indivq \
--ind_labels examples/capeverde_ind_labels.txt

# run with all plot types
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_output \
-f admixture \
--extension .indivq \
--plot_type all \
--ind_labels examples/capeverde_ind_labels.txt

# run with figures in .png format
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_output \
-f admixture \
--use_rep T \
--use_best_pair F \
--extension .indivq \
--plot_type graph \
--fig_format png \
--ind_labels examples/capeverde_ind_labels.txt

# run with customized mode detection
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_custom_output \
-f admixture \
--use_rep F \
--use_best_pair F \
--extension .indivq \
--cd_method custom \
--ind_labels examples/capeverde_ind_labels.txt

# run without reordering individuals within each label group, and plot unaligned modes as well
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_output_unordered_w_unaligned \
-f admixture \
--use_rep F \
--use_best_pair F \
--extension .indivq \
--reorder_within_group F \
--plot_type all \
--plot_unaligned T \
--ind_labels examples/capeverde_ind_labels.txt

# run without any reordering
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_output_unordered \
-f admixture \
--use_rep F \
--use_best_pair F \
--extension .indivq \
--regroup_ind F \
--reorder_within_group F \
--order_cls_by_label F \
--plot_type graph \
--ind_labels examples/capeverde_ind_labels.txt

# run without any individual reordering, but with cluster ordering by label
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_output_unordered_cls_order \
-f admixture \
--use_rep F \
--use_best_pair F \
--extension .indivq \
--regroup_ind F \
--reorder_within_group F \
--order_cls_by_label T \
--plot_type graph \
--ind_labels examples/capeverde_ind_labels.txt


