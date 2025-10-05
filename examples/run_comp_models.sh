mkdir -p examples/comp_models
model1=louvain # should not contain underscore or dash in the model name
model1_suffix=avg
model2=markov # should not contain underscore or dash in the model name
model2_suffix=avg
ls examples/capeverde_avg_output/modes_aligned/*_${model1_suffix}.Q > examples/comp_models/${model1}.qfilelist
for f in examples/capeverde_avg_output/modes_aligned/*_${model1_suffix}.Q; do [ -f "$f" ] && basename "$f" | sed "s/\_${model1_suffix}.Q$//" >> examples/comp_models/${model1}.qnamelist; done
ls examples/capeverde_custom_output/modes_aligned/*_${model2_suffix}.Q > examples/comp_models/${model2}.qfilelist
for f in examples/capeverde_custom_output/modes_aligned/*_${model2_suffix}.Q; do [ -f "$f" ] && basename "$f" | sed "s/\_${model2_suffix}.Q$//" >> examples/comp_models/${model2}.qnamelist; done
cp examples/capeverde_avg_output/modes/mode_stats.txt examples/comp_models/${model1}_mode_stats.txt
cp examples/capeverde_custom_output/modes/mode_stats.txt examples/comp_models/${model2}_mode_stats.txt

python -m clumppling.compModels \
--models ${model1} ${model2} \
--qfilelists examples/comp_models/${model1}.qfilelist examples/comp_models/${model2}.qfilelist \
--qnamelists examples/comp_models/${model1}.qnamelist examples/comp_models/${model2}.qnamelist \
--mode_stats_files examples/comp_models/${model1}_mode_stats.txt examples/comp_models/${model2}_mode_stats.txt \
--output examples/comp_models/capeverde_comp_${model1}_vs_${model2} 