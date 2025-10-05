# !/bin/bash

# default run with louvain mode detection, representative Q matrices for modes, and alignment across-K using the best mode pair
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_output \
-f admixture \
--extension .indivq \
--ind_labels examples/capeverde_ind_labels.txt

# run with average Q matrices for modes, and alignment across-K using the major mode pair
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_avg_output \
-f admixture \
--use_rep F \
--use_best_pair F \
--extension .indivq \
--ind_labels examples/capeverde_ind_labels.txt

# run with average Q matrices for modes, alignment across-K using the major mode pair, and all plot types
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_output \
-f admixture \
--use_rep F \
--use_best_pair F \
--extension .indivq \
--plot_type all \
--ind_labels examples/capeverde_ind_labels.txt

# run with average Q matrices for modes, alignment across-K using the major mode pair, and graph plot in .png format
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_output \
-f admixture \
--use_rep F \
--use_best_pair F \
--extension .indivq \
--plot_type graph \
--fig_format png \
--ind_labels examples/capeverde_ind_labels.txt

# run with customized mode detection, average Q matrices for modes, and alignment across-K using the major mode pair
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_custom_output \
-f admixture \
--use_rep F \
--use_best_pair F \
--extension .indivq \
--cd_method custom \
--ind_labels examples/capeverde_ind_labels.txt


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




