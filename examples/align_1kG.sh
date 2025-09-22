# 1kG from Pong
input_dir=examples/1kG-p3_clumppling_input
pop_label_file=$input_dir/population_labels.txt
output_dir=examples/1kG-p3_clumppling_output
mkdir -p $output_dir

# reproduce Pong's results
python -m clumppling -i $input_dir -o $output_dir \
-f admixture --extension .Q --ind_labels $pop_label_file \
--comm_min 1e-6 --test_comm F --cd_res 1.05 \
--plot_type all --regroup_ind T --alt_color T
