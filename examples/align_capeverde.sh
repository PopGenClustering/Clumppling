# !/bin/bash

# default run with louvain mode detection, representative Q matrices for modes, and alignment across-K using the best mode pair
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_output \
-f admixture \
--extension .indivq \
--ind_labels examples/capeverde_ind_labels.txt

# run with average Q matrices, average Q matrices for modes, and alignment across-K using the major mode pair
python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_avg_output \
-f admixture \
--use_rep F \
--use_best_pair F \
--extension .indivq \
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

