# !/bin/bash

# run with structure output files, all visualizations, and custom color map
python -m clumppling \
-i examples/chicken \
-o examples/chicken_output \
-f structure \
--extension _f \
--plot_type all \
--fig_format png \
--custom_cmap examples/custom_colors.txt 