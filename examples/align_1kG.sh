# 1kG from Pong
input_dir=examples/1kG-p3_clumppling_input
pop_label_file=$input_dir/population_labels.txt
uniq_pop_label_file=$input_dir/ordered_unique_population_labels.txt
output_dir=examples/1kG-p3_clumppling_output
mkdir -p $output_dir

# reproduce Pong's results
python -m clumppling -i $input_dir -o $output_dir \
-f admixture --extension .Q \
--ind_labels $pop_label_file --ordered_uniq_labels $uniq_pop_label_file \
--test_comm F --cd_res 1.05 --fig_format png \
--plot_type all --regroup_ind T --alt_color T
