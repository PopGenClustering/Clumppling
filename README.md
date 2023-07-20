# *Clumppling*

This is the GitHub repository for the program ***Clumppling*** (CLUster Matching and Permutation Program that uses integer Linear ProgrammING), a framework for aligning clustering results of population structure analysis.

Current version **v 0.1.0 (beta version)** 

> The current version has been tested on Windows 10 with Python 3.8 to 3.10, Ubuntu 20.04 LTS with Python 3.10, MacOS 13 with Python 3.11, and RHEL (Red Hat Enterprise Linux) 8.6 with Python 3.11.
>
> **Detailed instructions can be found in [the pdf Manual](Clumppling_Manual.pdf).**
>
> Questions and feedback are welcome. Contact Xiran Liu at xiranliu@stanford.edu.

**There are two ways to run *Clumppling*.**
1. You can run it **remotely** on the server, which does not require downloading or installing the program locally. The remote version provides the core functionalities of the program. Check out the [Remote Version](#Remote-Version) section.
2. You can download and install the Python package onto your local machine and run the program **locally**. The local version provides an extended list of functionalities (see [the pdf Manual](Clumppling_Manual.pdf) for details). Check out the [Local Version](#Local-Version) section.

---
# Remote Version
The remote version is available through an **online Colaboratory notebook**, which is a Jupyter notebook that runs in the cloud served by Google. If you are interested, more details about Colab notebooks can be found here: https://colab.google/.

There is **no need to download and install the program locally**.

To run *Clumppling* remotely, click on [this link](https://colab.research.google.com/drive/1PiM5pUKm9cx-dCz0YLWwaJcNcTQHyUm8#offline=true&sandboxMode=true) which will bring you to the notebook. \
Next, follow the instructions in the notebook. \
Click the run (little round-shaped buttons with a triangle in the middle) buttons next to each block on the left one by one. \
Upload input files (e.g., the example files provided [here](input)) as a zip folder, specify the input data format, and change input parameters (if needed) following the instructions. \
You will be able to download a zipped file containing the alignment results at the end of the notebook.

---
# Local Version
The local version requires downloading and installing the program to your local machine. 

## How to Install

### 1. Open a **command line interpreter** (i.e., a shell)
   * Linux and macOS users can use the built-in Terminal. 
   * For Windows users, you need to get a terminal. For example: 
     - After you follow Step 3 to install Conda, you can use the built-in *Anaconda Prompt* available from the Anaconda Navigator.
     - You may also use the built-in (for Windows 10) [*Windows PowerShell*](https://learn.microsoft.com/en-us/windows-server/administration/windows-commands/powershell). 
     - Or, you can use *Git Bash* after you install Git by downloading and running the executable Git installer from https://git-scm.com/download/win.

### 2. Install **Python** (Version 3.8 and above recommended)
   You can download the Python installer from https://www.python.org/downloads/.
   * For **Windows** users, go to https://www.python.org/downloads/windows/ to download the installer corresponding to your operating system, e.g., Windows installer (64-bit). Run the executable installer and check the box 'Add Python to environment variables' during the installation.
   * For **macOS** users, go to https://www.python.org/downloads/macos/ to download the macOS 64-bit universal2 installer and double-click on the *python-<version>-macosx.pkg* file to start the Python installer.
   * For **Linux** users, if Python is not pre-installed, you can install it via command lines (``sudo yum install -y python3`` for CentOS and Red Hat Linux and ``sudo apt-get install python3`` for all other Linux systems). 
   
   You can verify the installation by running 
   ````
   python --version
   ```` 
   in the command line interpreter, which should give you the version of the installed Python.

### 3. Install conda and create a virtual environment 
   Go to https://www.anaconda.com/download to download the conda installer and run the installer. Conda is a popular package management system and environment management system.
   
   > A virtual environment is a Python environment such that the Python interpreter, libraries and scripts installed into it are isolated from those installed in other virtual environments, and (by default) any libraries installed in a “system” Python, i.e., one which is installed as part of your operating system"
   
   Using a virtual environment helps to keep the dependencies required by different projects separate and avoid conflicts between projects. 
   
   Create a virtual environment named ``clumppling-env`` by type the following command in the command-line interpreter
   ````
   conda create -n clumppling-env python
   ````
   Activate the virtual environment by
   ````
   conda activate clumppling-env 
   ````
### 4. Install the *Clumppling* package 
   **(1) Install the package**  

   If you have [Git](https://git-scm.com/) installed, run the following command:
   ````
   pip install git+https://github.com/PopGenClustering/Clumppling
   ````
   If you don't have Git, run the following command:
   ````
   pip install https://github.com/PopGenClustering/Clumppling/archive/master.zip
   ````

   **(2) Download the example files from [the input directory](input) in the GitHub repository** \
   For each example dataset, unzip the files into a folder with the same name, and put it under a folder called "input" under a path of your choice. 
   
   For example, you may put it under a path ``C:/Users/YOUR_USER_NAME/Clumppling``. Then the path for the example files should be ``C:/Users/YOUR_USER_NAME/Clumppling/input``, and the path for the Cape Verde example data should be ``C:/Users/YOUR_USER_NAME/Clumppling/input/capeverde``.
   
   More will be discussed in the section *How to Run (with an example)*.
      
### 5. Check whether the installation is successful
   Run the following command:
   ````
   python -m clumppling -h
   ````
   If the installation was successful, you should see the usage of the program in the command window which tells you the required and optional arguments to the program. It should look like:
   ````
    usage: __main__.py [-h] -i INPUT_PATH -o OUTPUT_PATH -f INPUT_FORMAT [-v VIS] [--cd_param CD_PARAM]
                       [--use_rep USE_REP] [--merge_cls MERGE_CLS] [--cd_default CD_DEFAULT] [--plot_modes PLOT_MODES]
                       [--plot_modes_withinK PLOT_MODES_WITHINK] [--plot_major_modes PLOT_MAJOR_MODES]
                       [--plot_all_modes PLOT_ALL_MODES] [--custom_cmap CUSTOM_CMAP]
    
    required arguments:
      -i INPUT_PATH, --input_path INPUT_PATH
                            path to load input files
      -o OUTPUT_PATH, --output_path OUTPUT_PATH
                            path to save output files
      -f INPUT_FORMAT, --input_format INPUT_FORMAT
                            input data format
    
    optional arguments:
      -v VIS, --vis VIS     whether to generate visualization: 0 for no, 1 for yes (default)
      --cd_param CD_PARAM   the parameter for community detection method (default 1.0)
      --use_rep USE_REP     whether to use representative replicate as mode consensus: 0 for no (default), 1 for yes
      --merge_cls MERGE_CLS
                            whether to merge all pairs of clusters to align K+1 and K: 0 for no (default), 1 for yes
      --cd_default CD_DEFAULT
                            whether to use default community detection method (Louvain): 0 for no, 1 for yes (default)
      --plot_modes PLOT_MODES
                            whether to display aligned modes in structure plots over a multipartite graph: 0 for no, 1 for
                            yes (default)
      --plot_modes_withinK PLOT_MODES_WITHINK
                            whether to display modes for each K in structure plots: 0 for no (default), 1 for yes
      --plot_major_modes PLOT_MAJOR_MODES
                            whether to display all major modes in a series of structure plots: 0 for no (default), 1 for
                            yes
      --plot_all_modes PLOT_ALL_MODES
                            whether to display all aligned modes in a series of structure plots: 0 for no (default), 1 for
                            yes
      --custom_cmap CUSTOM_CMAP
                            customized colormap as a comma-separated string of hex codes for colors: if empty (default),
                            using the default colormap, otherwise use the user-specified colormap
   ````

## Input Arguments
The main module takes in three required arguments and several optional ones. The required arguments are
* ``-i`` (``--input_path``) path to load input files
* ``-o`` (``--output_path``) path to save output files
* ``-f`` (``--input_format``) input data format. This choice must be one of "structure", "admixture", "fastStructure", and "generalQ".

The optional arguments are ``cd_param``, ``use_rep``, ``merge_cls``, ``cd_default``, ``plot_modes``, ``plot_modes_withinK``, ``plot_major_modes``, ``plot_all_modes``, and ``custom_cmap``. **Detailed explanation of these arguments can be found in the [Manual](Clumppling_Manual.pdf).**

    
## How to Run (with an example)
As a quick start, let's use the [Cape Verde data](https://doi.org/10.1016/j.cub.2017.07.002) as an example. The data files are available in the zip file [input/capeverde.zip](input/capeverde.zip). 

1. **Ensure that the data files have been successfully downloaded and put under the right directory.** Following Step 4(2) from above, you should have already downloaded the zip file "capeverde.zip" and extracted the data files inside the zip file to a folder named "capeverde" under a directory called "input" that you created. \
For example, the data files are extracted to ``C:/Users/YOUR_USER_NAME/Clumppling/input/capeverde``.

2. **Ensure that the current path is the parent directory of the "input" folder.** In your command-line interpreter, make sure that you navigate to the directory where the folder "input" locates. For the above example, you should be in the directory ``C:/Users/YOUR_USER_NAME/Clumppling``. 
   
   If you do not know your current path, you can run ``pwd`` to see it in the command window. 
   
   The Cape Verde data files should locate in ``input/capeverde`` under your current path. I.e., When you run ``ls input/capeverde``, you should see a list of Cape Verde data files *(CAPEVERDE_Rep1.2.indivq, CAPEVERDE_Rep1.3.indivq, etc.)*

3. **Run the program** under the default setting
   ````
   python -m clumppling -i input/capeverde -o output/capeverde -f admixture 
   ````
   The outputs will be saved in ``output/capeverde`` under your current and a zipped file of the same name will also be generated. For example, if you are in the directory ``C:/Users/YOUR_USER_NAME/Clumppling``, then the input will be saved in ``C:/Users/YOUR_USER_NAME/Clumppling/output/capeverde`` and zipped in ``C:/Users/YOUR_USER_NAME/Clumppling/output/capeverde.zip``.
   
   ***Notes:***
   * Make sure that the output path (e.g. ''output/capeverde'') is different from the input path (e.g. ''input/capeverde'')! Otherwise, the input directory will be overwritten by the outputs and this causes issues.
   * If your data files are in a directory different than the given example, make sure to specify the corresponding input path after ``-i``.
   * You can also change the parameters from those in the default setting. Here is another example for the chicken example (data available in the zip file [input/chicken.zip](input/chicken.zip)) where we specify the parameters ``cd_param`` and ``custom_cmap``:
      ````
      python -m clumppling -i input/chicken -o output/chicken_color -f structure --cd_param 1.05 --custom_cmap #D65859,#00AAC1,#01C0F6,#FDF0C4,#F1B38C,#AAD6BD,#6BB582,#B5DDF7,#AE8557,#FCEC73,#A4A569,#4264AC,#A1CDB2,#DE9D5D,#D9439A,#ABB2BA,#8775B3,#B3865C,#DADDE6,#E7BDD1,#FF9999
      ````
      Here we are setting the community detection parameter ``cd_param`` to be 1.05. This parameter controls the resolution, i.e., the size of the detected communities, of the Louvain algorithm used for community detection. The default value is 1.0, and a larger value will result in smaller communities and more of them, thus detecting more modes. If the default parameter value is giving undesirable mode detection results, for instance, assigning each replicate to its own modes or assigning all obviously different replicates to a single mode, you may vary this parameter to get a desired size and number of the modes. \
      We are also providing a custom colormap (``custom_cmap``) for visualization of the alignment results. The clusters will be colored based on colors specified by the list of comma-separated HEX codes.

   * More example commands to call the program can be found under [the scripts directory](scripts) in the GitHub repository.


**Detailed instructions to run the program can be found in the [Manual](Clumppling_Manual.pdf).**

## References
*Liu, Kopelman, & Rosenberg (2023). Clumppling: cluster matching and permutation program with integer linear programming. \[in prep\]*

The Cape Verde data used as the example comes from: \
*Verdu, P., Jewett, E. M., Pemberton, T. J., Rosenberg, N. A., & Baptista, M. (2017). Parallel trajectories of genetic and linguistic admixture in a genetically admixed creole population. Current Biology, 27(16), 2529-2535. [https://doi.org/10.1016/j.cub.2017.07.002.](https://www.sciencedirect.com/science/article/pii/S096098221730859X)*

The chicken data used as the example comes from: \
*Rosenberg, N. A., Burke, T., Elo, K., Feldman, M. W., Freidlin, P. J., Groenen, M. A., ... & Weigend, S. (2001). Empirical evaluation of genetic clustering methods using multilocus genotypes from 20 chicken breeds. Genetics, 159(2), 699-713. [https://doi.org/10.1093/genetics/159.2.699](https://doi.org/10.1093/genetics/159.2.699)*

## Acknowledgements
We thank Egor Lappo for helping with the packaging of the program. 
We thank Egor Lappo, Daniel Cotter, Maike Morrison, Chloe Shiff, and Juan Esteban Rodriguez Rodriguez for helping with the testing of the program.
