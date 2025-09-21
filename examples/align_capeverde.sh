python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_output \
-f admixture \
--extension .indivq \
--ind_labels examples/capeverde_ind_labels.txt

python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_avg_output \
-f admixture \
--use_rep F \
--use_best_pair F \
--extension .indivq \
--ind_labels examples/capeverde_ind_labels.txt

python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_custom_output \
-f admixture \
--use_rep F \
--use_best_pair F \
--extension .indivq \
--cd_method custom \
--ind_labels examples/capeverde_ind_labels.txt

