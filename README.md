# Clumppling

This is the GitHub repository for the program ***Clumppling*** (CLUster Matching and Permutation Program that uses integer Linear programmING), a framework for aligning mixed-membership clustering results of population structure analysis.

Current version **v 1.0.1** (Last update: Aug 14, 2025)

## Features

- **Flexible input parsing**
- **Cluster alignment within and across K values**
- **Mode detection in clustering results**
- **Visualization of alignment patterns and aligned modes**
- **Modular design for easy integration**

## Installation

```bash
pip install clumppling
```

Run the following command:
   ````bash
   python -m clumppling -h
   ````
If the installation was successful, you should see the usage of the program in the command window. The usage tells you the required and optional arguments to the program. It should look like:
````bash
usage: __main__.py [-h] -i INPUT -o OUTPUT -f {generalQ,admixture,structure,fastStructure} [-v VIS] [--custom_cmap CUSTOM_CMAP] [--plot_type {graph,list,withinK,major,all}]
                [--include_cost INCLUDE_COST] [--include_label INCLUDE_LABEL] [--ind_labels IND_LABELS] [--reorder_ind REORDER_IND] [--reorder_by_max_k REORDER_BY_MAX_K]
                [--order_cls_by_label ORDER_CLS_BY_LABEL] [--extension EXTENSION] [--skip_rows SKIP_ROWS] [--remove_missing REMOVE_MISSING]
                [--cd_method {louvain,leiden,infomap,markov_clustering,label_propagation,walktrap,custom}] [--cd_res CD_RES] [--test_comm TEST_COMM] [--comm_min COMM_MIN]
                [--comm_max COMM_MAX] [--merge MERGE] [--use_rep USE_REP] [--use_best_pair USE_BEST_PAIR]

Clumppling: a tool for cluster matching and permutation program with integer linear programming

required arguments:
-i INPUT, --input INPUT
                        Input file path
-o OUTPUT, --output OUTPUT
                        Output file directory
-f {generalQ,admixture,structure,fastStructure}, --format {generalQ,admixture,structure,fastStructure}
                        File format

optional arguments:
-v VIS, --vis VIS     Whether to generate figure(s): True (default)/False
--custom_cmap CUSTOM_CMAP
                        Customized colormap as a comma-separated string of hex codes for colors: if empty (default), using the default colormap, otherwise use the user-specified
                        colormap
--plot_type {graph,list,withinK,major,all}
                        Type of plot to generate: 'graph' (default), 'list', 'withinK', 'major', 'all'
--include_cost INCLUDE_COST
                        Whether to include cost values in the graph plot: True (default)/False
--include_label INCLUDE_LABEL
                        Whether to include individual labels in the plot: True (default)/False
--ind_labels IND_LABELS
                        A plain text file containing individual labels (one per line) (default: last column from labels in input file, which consists of columns [0, 1, 3]
                        separated by delimiter)
--reorder_ind REORDER_IND
                        Whether to reorder individuals within each label group in the plot: True (default)/False
--reorder_by_max_k REORDER_BY_MAX_K
                        Whether to reorder individuals based on the major mode with largest K: True (default)/False (based on the major mode with smallest K)
--order_cls_by_label ORDER_CLS_BY_LABEL
                        Whether to reorder clusters based on total memberships within each label group in the plot: True (default)/False (by overall total memberships)
--extension EXTENSION
                        Extension of input files
--skip_rows SKIP_ROWS
                        Skip top rows in input files
--remove_missing REMOVE_MISSING
                        Remove individuals with missing data: True (default)/False
--cd_method {louvain,leiden,infomap,markov_clustering,label_propagation,walktrap,custom}
                        Community detection method to use (default: louvain)
--cd_res CD_RES       Resolution parameter for the default Louvain community detection (default: 1.0)
--test_comm TEST_COMM
                        Whether to test community structure (default: True)
--comm_min COMM_MIN   Minimum threshold for cost matrix (default: 1e-4)
--comm_max COMM_MAX   Maximum threshold for cost matrix (default: 1e-2)
--merge MERGE         Whther to merge two clusters when aligning K+1 to K (default: True)
--use_rep USE_REP     Use representative modes (alternative: average): True (default)/False
--use_best_pair USE_BEST_PAIR
                        Use best pair as anchor for across-K alignment (alternative: major): True (default)/False
````

## Usage
Examples:
````bash
python -m clumppling.clumppling \
        -i tests/test1/input \
        -o tests/test1/output \
        -f generalQ 
````

## Main Function

**Detailed explanation will be available in the [Manual](Clumppling_Manual.pdf).**

### Input Arguments
The main module takes in three required arguments and several optional ones. The required arguments are
* ``-i`` (``--input``) path to load input files
* ``-o`` (``--output``) path to save output files
* ``-f`` (``--format``) input data format. This choice must be one of "generalQ", "admixture", "structure", or "fastStructure".

The optional arguments are 
* for input parsing: ``extension``, ``skip_rows``, ``remove_missing``
* for community detection: ``cd_method``,  ``cd_res``, ``test_comm``, ``comm_min``, ``comm_max``
* for alignment across-K: ``merge``, ``use_rep``,``use_best_pair``
* for figure generation: ``vis``, ``plot_type``,``include_cost``, ``include_label``, ``ind_labels``, ``custom_cmap``, ``reorder_ind``, ``reorder_by_max_k``, ``order_cls_by_label``.

### How to Run (with an example)
As a quick start, let's use the [Cape Verde data](https://doi.org/10.1016/j.cub.2017.07.002) and the [chicken data](https://doi.org/10.1093/genetics/159.2.699) as examples. The data files are available in the zip files [examples/capeverde.zip](examples/capeverde.zip) and [examples/chicken.zip](examples/chicken.zip). 

The *.indivq* files for Cape Verde data contains the column indicating their population indices. The corresponding population labels are provided separately in the file [examples/capeverde_ind_labels.txt](examples/capeverde_ind_labels.txt).

A file with custom colors is also provided at [examples/custom_colors.txt](examples/custom_colors.txt) for use in examples.

1. **Ensure that the data files have been successfully downloaded and put under the right directory.** 
Download the example files from [the directory "examples"](examples) in the GitHub repository. For each example dataset, unzip the files into the folder with the same name as the zip file. 
   
2. **Ensure that the current path is the correct directory.** By default, you should be in the parent directory of the "examples" folder, i.e., in your command-line interpreter, make sure that you navigate to the directory where the folder "exmaples" is located. Alternatively, update the paths correspondingly in the following example scripts.

3. **Run the program** on the Cape Verde data under the default setting, with user-provided individual labels:

    ````bash
    python -m clumppling \
    -i examples/capeverde \
    -o examples/capeverde_output \
    -f admixture \
    --extension .indivq \
    --ind_labels examples/capeverde_ind_labels.txt
    ````

    The outputs will be saved in "examples/capeverde_output" under your current directory and a zipped file of the same name will also be generated and zipped in ``examples/capeverde_output.zip``.

    Similarly, you can run the program on the chicken data as follows:

    ````bash
    python -m clumppling \
    -i examples/chicken \
    -o examples/chicken_output \
    -f structure \
    --extension _f \
    --plot_type all \
    --use_best_pair 0 \
    --custom_cmap examples/custom_colors.txt 
    ````

    The outputs will be saved in "examples/chicken_output" under your current directory and a zipped file of the same name will also be generated and zipped in ``examples/chicken_output.zip``.
    
These commands are also provided in the [example script for running Clumppling on Cape Verde data (Admixture .indviq files)](examples/align_capeverde.sh) and [example script for running Clumppling on chicken data (Structure _f files)](examples/align_chicken.sh). 


### Outputs
The output folder will contain the following structure (see `examples/capeverde_output` for reference after finishing running the example; suppose `use_rep=True`):

```
output/
├── input/
│   ├── input_meta.txt
│   ├── input_labels.txt (optional)
│   ├── 1_K2R1.Q
│   └── ...
├── alignment_withinK/
│   ├── K2.txt
│   ├── K3.txt
│   └── ...
├── modes/
│   ├── mode_alignment.txt
│   ├── mode_stats.txt
│   ├── mode_average_stats.txt
│   ├── K2M1_avg.Q
│   ├── K2M1_rep.Q
│   └── ...
├── modes_aligned/
│   ├── all_modes_alignment_rep.txt
│   ├── K2M1_rep.Q
│   └── ...
├── alignment_acrossK/
│   ├── alignment_acrossK_rep.txt
│   ├── best_pairs_acrossK_rep.txt
│   └── major_pairs_acrossK_rep.txt
├── visualization/
│   ├── colorbar.png
│   ├── alignment_pattern_rep.png
│   ├── all_modes_graph_rep.png
│   └── ...
└── clumppling.log
```

- `aligned_modes/`: Contains files with clusters aligned within each K.
- `modes/`: Contains detected modes for each K.
- `acrossK_alignment/`: Contains results of cluster alignment across different K values.
- `plots/`: Contains generated visualizations (structure plots, alignment graphs).
- `logs/`: Contains log files from the run.

File names and subfolders may vary depending on your input and options.

## Submodules

Each submodule is callable independently.

### `parseInput`

Handles reading and parsing input files containing clustering results. Supports various formats and prepares data for downstream analysis.

Example:
````bash
python -m clumppling.parseInput \
-i examples/submodules/input \
-o examples/submodules/output \
-f generalQ 
````

### `aligneWithinK`

Aligns clusters within a single value of K to ensure consistent labeling and facilitate comparison across replicates.

Example:
````bash
python -m clumppling.alignWithinK \
--qfilelist examples/submodules/K3.qfilelist \
-o examples/submodules/K3_aligned.txt

python -m clumppling.alignWithinK \
--qfilelist examples/submodules/K5.qfilelist \
-o examples/submodules/K5_aligned.txt
````

### `detectMode`

Detects modes (distinct clustering solutions) among multiple runs for a given K, helping to identify stable and alternative solutions.

Example:
````bash
python -m clumppling.detectMode \
--align_res examples/submodules/K3_aligned.txt \
--qfilelist examples/submodules/K3.qfilelist \
--qnamelist examples/submodules/K3.qnamelist \
--cd_method markov_clustering \
-o examples/submodules/K3_modes

python -m clumppling.detectMode \
--align_res examples/submodules/K5_aligned.txt \
--qfilelist examples/submodules/K5.qfilelist \
--qnamelist examples/submodules/K5.qnamelist \
--cd_method markov_clustering \
-o examples/submodules/K5_modes
````

### `alignAcrossK`

Aligns clusters across different values of K, enabling tracking of cluster membership changes as K varies.

Example:
````bash
# prepare mode files
for K in 3 5; do
    SRC=examples/submodules/K${K}_modes
    DST=examples/submodules/K3K5_modes
    mkdir -p $DST
    for f in "$SRC"/*.Q; do
        if [ -f "$f" ]; then
            cp "$f" "$DST/K${K}$(basename "$f")"
        fi
    done
done
# generate list of mode files and names
ls examples/submodules/K3K5_modes/*_rep.Q > examples/submodules/K3K5_modes/K3K5_modes.qfilelist
for f in examples/submodules/K3K5_modes/*_rep.Q; do [ -f "$f" ] && basename "$f" | sed 's/\_rep.Q$//' >> examples/submodules/K3K5_modes/K3K5_modes.qnamelist; done

# run clumppling
python -m clumppling.alignAcrossK \
--qfilelist examples/submodules/K3K5_modes/K3K5_modes.qfilelist \
--qnamelist examples/submodules/K3K5_modes/K3K5_modes.qnamelist \
-o examples/submodules/K3K5_acrossK_output
````

### Visualizations
For generating figures, see [examples/plot_submodules.py](examples/plot_submodules.py) as an example. Run
````bash
python examples/plot_submodules.py 
````

All commands are also provided in [examples/run_submodules.sh](examples/run_submodules.sh). 


### `compModels`

Compare clustering results from different models, with potentially different K values for the results of each model.

````bash
python -m clumppling.compModels \
--models ${model1} ${model2} \
--qfilelists examples/comp_models/${model1}.qfilelist examples/comp_models/${model2}.qfilelist \
--qnamelists examples/comp_models/${model1}.qnamelist examples/comp_models/${model2}.qnamelist \
--output examples/comp_models/capeverde_comp_${model1}_vs_${model2} 
````

Example: see [examples/run_comp_models.sh](examples/run_comp_models.sh). 


## License

MIT License

## References
*Liu, X., Kopelman, N. M., & Rosenberg, N. A. (2024). Clumppling: cluster matching and permutation program with integer linear programming. Bioinformatics, 40(1), btad751. [https://doi.org/10.1093/bioinformatics/btad751](https://doi.org/10.1093/bioinformatics/btad751)*

The Cape Verde data used as the example comes from: \
*Verdu, P., Jewett, E. M., Pemberton, T. J., Rosenberg, N. A., & Baptista, M. (2017). Parallel trajectories of genetic and linguistic admixture in a genetically admixed creole population. Current Biology, 27(16), 2529-2535. [https://doi.org/10.1016/j.cub.2017.07.002.](https://www.sciencedirect.com/science/article/pii/S096098221730859X)*

The chicken data used as the example comes from: \
*Rosenberg, N. A., Burke, T., Elo, K., Feldman, M. W., Freidlin, P. J., Groenen, M. A., ... & Weigend, S. (2001). Empirical evaluation of genetic clustering methods using multilocus genotypes from 20 chicken breeds. Genetics, 159(2), 699-713. [https://doi.org/10.1093/genetics/159.2.699.](https://doi.org/10.1093/genetics/159.2.699)*

## Acknowledgements
We thank Egor Lappo for helping with the packaging of the program. 
We thank Egor Lappo, Daniel Cotter, Maike Morrison, Chloe Shiff, and Juan Esteban Rodriguez Rodriguez for helping with the testing of the program.

## Version Update History
Version 0.3.2: 
- Fix the bug in plotting when all replicates have the same K.
- Add a check (merged from branch) to exclude loading K=1 replicates.

Version 1.0.0 (major updates): 
- Modularize each step.
- Add visualization of alignment patterns.
- Add input parsing features:
    * use `extension` to specify the file extension of the input files.
    * use `skip_rows` to specify number of rows to skip from the input files.
    * use `remove_missing` to choose whether to remove individuals with missing data (clusters).
- Add more flexibility in algorithmic settings: 
    * `test_comm`: whether to test for community structure during mode detection, as well as extreme values for determining if nodes fall into communities (`comm_min` and `comm_max`). 
    * `use_best_pair`: whether to align across-K using the best pair of modes or the pair of major modes as the anchor.
    * keep the features `merge` and `use_rep`.
- Use the package `cdlib` for community detection. Change `cd_method` to multiple choices (default: 'louvain') and move `cd_custom` as the choice 'custom'. 
- Add more flexibility in plotting settings:
    * `plot_type`: which plot(s) to generate: 'all', 'graph' (default), 'list', 'major', or 'withinK'.
    * `include_cost`: include edges indicating alignment costs in the graph of structure plots.
    * `include_label`: whether to include group labels of individuals (if available) on the x-axis and draw corresponding vertical lines in the structure plots separating groups.
    * `ind_labels`: accept user-specified individual labels from a file.

> Questions and feedback are welcome.
> Contact the author at ``xiranliu at stanford dot edu``.