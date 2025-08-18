mkdir -p examples/comp_models
model1=louvain # should not contain dash in the model name
model1_suffix=avg
model2=markov # should not contain dash in the model name
model2_suffix=avg
ls examples/capeverde_output/modes_aligned/*_${model1_suffix}.Q > examples/comp_models/${model1}.qfilelist
for f in examples/capeverde_output/modes_aligned/*_${model1_suffix}.Q; do [ -f "$f" ] && basename "$f" | sed "s/\_${model1_suffix}.Q$//" >> examples/comp_models/${model1}.qnamelist; done
ls examples/capeverde_custom_output/modes_aligned/*_${model2_suffix}.Q > examples/comp_models/${model2}.qfilelist
for f in examples/capeverde_custom_output/modes_aligned/*_${model2_suffix}.Q; do [ -f "$f" ] && basename "$f" | sed "s/\_${model2_suffix}.Q$//" >> examples/comp_models/${model2}.qnamelist; done

python -m clumppling.compModels \
--models ${model1} ${model2} \
--qfilelists examples/comp_models/${model1}.qfilelist examples/comp_models/${model2}.qfilelist \
--qnamelists examples/comp_models/${model1}.qnamelist examples/comp_models/${model2}.qnamelist \
--output examples/comp_models/capeverde_comp_${model1}_vs_${model2} 