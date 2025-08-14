# Clumppling

This is the GitHub repository for the program ***Clumppling*** (CLUster Matching and Permutation Program that uses integer Linear programmING), a framework for aligning mixed-membership clustering results of population structure analysis.

Current version **v 1.0.0** (Last update: Aug 13, 2025)

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
   usage: __main__.py [-h] -i INPUT -o OUTPUT -f {generalQ,admixture,structure,fastStructure} [-v VIS] [--custom_cmap CUSTOM_CMAP] [--plot_type {graph,list,withinK,major,all}] [--include_cost INCLUDE_COST]
                   [--include_label INCLUDE_LABEL] [--ind_labels IND_LABELS] [--extension EXTENSION] [--skip_rows SKIP_ROWS] [--remove_missing REMOVE_MISSING]
                   [--cd_method {louvain,leiden,infomap,markov_clustering,label_propagation,walktrap,custom}] [--cd_res CD_RES] [--test_comm TEST_COMM] [--comm_min COMM_MIN] [--comm_max COMM_MAX]
                   [--merge MERGE] [--use_rep USE_REP] [--use_best_pair USE_BEST_PAIR]

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
                            Customized colormap as a comma-separated string of hex codes for colors: if empty (default), using the default colormap, otherwise use the user-specified colormap
    --plot_type {graph,list,withinK,major,all}
                            Type of plot to generate: 'graph' (default), 'list', 'withinK', 'major', 'all'
    --include_cost INCLUDE_COST
                            Whether to include cost values in the graph plot: True (default)/False
    --include_label INCLUDE_LABEL
                            Whether to include individual labels in the plot: True (default)/False
    --ind_labels IND_LABELS
                            A plain text file containing individual labels (one per line) (default: last column from labels in input file, which consists of columns [0, 1, 3] separated by delimiter)
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
* for figure generation: ``vis``, ``plot_type``,``include_cost``, ``include_label``, ``ind_labels``, ``custom_cmap``.

### How to Run (with an example)

### Outputs
The output folder will contain the following structure (see `tests/test1/output` for reference; suppose `use_rep=True`):

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
-i tests/test1/input \
-o tests/test1/output \
-f generalQ 
````

### `aligneWithinK`

Aligns clusters within a single value of K to ensure consistent labeling and facilitate comparison across replicates.

Example:
````bash
python -m clumppling.alignWithinK \
--qfilelist tests/test1/K5.qfilelist \
-o tests/test1/K5_aligned.txt
````

### `detectMode`

Detects modes (distinct clustering solutions) among multiple runs for a given K, helping to identify stable and alternative solutions.

Example:
````bash
python -m clumppling.detectMode \
--align_res tests/test1/K5_aligned.txt \
--qnamelist tests/test1/K5.qnamelist \
--cd_method markov_clustering \
-o tests/test1/K5_modes
````

### `alignAcrossK`

Aligns clusters across different values of K, enabling tracking of cluster membership changes as K varies.

Example:
````bash
python -m clumppling.alignAcrossK \
--qfilelist tests/test1/K3K5_modes/K3K5_modes.qfilelist \
--qnamelist tests/test1/K3K5_modes/K3K5_modes.qnamelist \
-o tests/test1/K3K5_modes/output
````

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