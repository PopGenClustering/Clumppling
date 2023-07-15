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
     - After you follow Step 3 to install Conda, you can use the built-in Anaconda Prompt available from the Anaconda Navigator.
     - Or, you can use Git Bash after you install Git by downloading and running the executable Git installer from https://git-scm.com/download/win.
     - You may also use the built-in (for Windows 10) [Windows PowerShell](https://learn.microsoft.com/en-us/windows-server/administration/windows-commands/powershell). 

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
   
   Create a virtual environment named ``clumppling-env`` by
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
   For each example dataset, unzip the files into a folder with the same name, and put it under a folder under your current path called "input". More will be discussed in the section *How to Run (with an example)*.
      
### 5. Check whether the installation is successful
   Run the following command:
   ````
   python -m clumppling -h
   ````
   If the installation was successful, you should see the usage of the program in the command window.

## Input Arguments
The main module takes in three required arguments and several optional ones. The required arguments are
* ``-i`` (``--input_path``) path to load input files
* ``-o`` (``--output_path``) path to save output files
* ``-f`` (``--input_format``) input data format. This choice must be one of "structure", "admixture", "fastStructure", and "generalQ".

The optional arguments are ``cd_param``, ``use_rep``, ``merge_cls``, ``cd_default``, ``plot_modes``, ``plot_modes_withinK``, ``plot_major_modes``, ``plot_all_modes``, and ``custom_cmap``. **Detailed explanation of these arguments can be found in the [Manual](Clumppling_Manual.pdf).**

    
## How to Run (with an example)
As a quick start, let's use the [Cape Verde data](https://doi.org/10.1016/j.cub.2017.07.002) as an example. The data files are available in the zip file [input/capeverde.zip](input/capeverde.zip). 

1. **Download the zip file** if you haven't yet (click the button on the top right to download raw file).

2. **Extract files** from the zip to a folder named "capeverde". Create a folder "input" under your current path and **put the "capeverde" folder in the "input" folder**.
   
   If you do not know your current path, you may run ``pwd`` to see it in the command window. 
   
   The chicken data files should be under ``input/capeverde``. I.e., When you run ``ls input/capeverde``, you should see a list of Cape Verde data files *(CAPEVERDE_Rep1.2.indivq, CAPEVERDE_Rep1.3.indivq, etc.)*

3. **Run the program** under the default setting
   ````
   python -m clumppling -i input/capeverde -o output/capeverde -f admixture 
   ````
   The outputs will be saved in the directory ``output/capeverde`` and a zipped file of the same name will also be generated.
   
   ***Notes:***
   * The output path (e.g. ''output/capeverde'') has to be different from the input path (e.g. ''input/capeverde'')! Otherwise, the input path will be overwritten by the outputs and this causes issues.
   * If your data files are in a different directory, make sure to specify the corresponding input path after ``-i``.
   * You can also change the parameters from those in the default setting. Here is another example for the chicken example (data available in the zip file [input/chicken.zip](input/chicken.zip)) where we specify the parameters ``cd_param`` and ``custom_cmap``:
      ````
      python -m clumppling -i input/chicken -o output/chicken_color -f structure --cd_param 1.05 --custom_cmap #D65859,#00AAC1,#01C0F6,#FDF0C4,#F1B38C,#AAD6BD,#6BB582,#B5DDF7,#AE8557,#FCEC73,#A4A569,#4264AC,#A1CDB2,#DE9D5D,#D9439A,#ABB2BA,#8775B3,#B3865C,#DADDE6,#E7BDD1,#FF9999
      ````
      Here we are setting the community detection parameter ``cd_param`` to be 1.05 and providing a custom colormap (``custom_cmap``) for visualization of the alignment results. 

   * More example commands to call the program can be found under [the scripts directory](scripts) in the GitHub repository.


**Detailed instructions to run the program can be found in the [Manual](Clumppling_Manual.pdf).**

## References
*Liu, Kopelman, & Rosenberg (2023). Clumppling: cluster matching and permutation program with integer linear programming. \[in prep\]*

The Cape Verde data used as the example comes from: \
*Verdu, P., Jewett, E. M., Pemberton, T. J., Rosenberg, N. A., & Baptista, M. (2017). Parallel trajectories of genetic and linguistic admixture in a genetically admixed creole population. Current Biology, 27(16), 2529-2535. [https://doi.org/10.1016/j.cub.2017.07.002.](https://www.sciencedirect.com/science/article/pii/S096098221730859X)*

## Acknowledgements
We thank Egor Lappo for helping with the packaging of the program. 
We thank Egor Lappo, Daniel Cotter, Maike Morrison, Chloe Shiff, and Juan Esteban Rodriguez Rodriguez for helping with the testing of the program.
