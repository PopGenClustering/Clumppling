python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_output \
-f admixture \
--extension .indivq \
--use_rep F \
--plot_type all \
--reorder_ind T \
--reorder_by_max_k T \
--order_cls_by_label T \
--ind_labels examples/capeverde_ind_labels.txt

python -m clumppling \
-i examples/capeverde \
-o examples/capeverde_custom_output \
-f admixture \
--extension .indivq \
--cd_method custom \
--use_rep F \
--ind_labels examples/capeverde_ind_labels.txt

