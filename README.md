# Clumppling

This is the GitHub repository for the program ***Clumppling*** (CLUster Matching and Permutation Program that uses integer Linear ProgrammING), a new framework for aligning clustering results of population structure analysis.

> Current version **v 0.1.0 (beta version)**.
> 
> The current version has been tested on Windows 10 with Python 3.8.0.

## How to Install
1. Install [Python3](https://www.python.org/downloads/). Versions above 3.7 are recommended. 
2. Install the *Clumppling* package. This can be done by either one of the following two ways:
    * Download the package files at https://github.com/PopGenClustering/Clumppling/archive/refs/heads/master.zip and unzip to a local directory. Remember to rename the folder named ``Clumppling-master`` to ``Clumppling``.
      
      Or, directly clone this repository to a local directory via ``git clone https://github.com/PopGenClustering/Clumppling.git``.
      
      Next, navigate into the directory ``Clumppling`` and run the following command
      ````
      pip install -e .
      ````

    * Install the package directly from GitHub via 
       ````
       pip install git+https://github.com/PopGenClustering/Clumppling
       ```` 
       However, **this will not get you the example data files automatically**. To run examples, you will need to download the example files from the ```input``` repository separately. 
3. Check the installation by running
      ````
      python clumppling -h
      ````
      If the installation was successful, you should be prompted by the usage of this function.
      
      Note that on some system, you may need to call ```python3``` instead of ```python```. You can easily set the alias by ```alias python=python3```.

## Input Arguments
The main module takes in three required arguments and several optional ones. The required arguments are
* ``-i`` (``--input_path``) path to load input files
* ``-o`` (``--output_path``) path to save output files
* ``-f`` (``--input_format``) input data format, have to be one of "structure", "admixture", "fastStructure", and "generalQ"

The optional arguments are ``cd_param``, ``use_rep``, ``merge_cls``, ``cd_default``, ``plot_modes``, ``plot_modes_withinK``, ``plot_major_modes``, ``plot_all_modes``, and ``custom_cmap``. **Detailed explanation of these arguments can be found in the Manual.**

    
## How to Run (Example: inferred population structure from chicken microsatellite data)
As a quick start, let's use the chicken data as an example. 

The data files should be available in the zip file ```input/chicken.zip```. 

* Unzip all data files in it to a folder named ```chicken``` under the same directory ```input```.

* Navigate to the directory ```Clumppling``` (where the ```README.md``` file is located) if you choose to download and install. Or, navigate to an arbitrary local directory.

* Once you are in the directory ```Clumppling``` and have downloaded and unzipped the chicken data files (in ``input/chicken``), run
   ````
   python clumppling -i input/chicken -o output/chicken -f structure 
   ````
   This will run the program with default parameters on the clustering results from chicken dataset. The outputs will be saved in the directory ``output/chicken`` and a zipped file with the same contents and a same name will also be created.

* If you are in a different directory, or if you put data files in a different directory, just make sure to specify the corresponding input path after ``-i``. 

* If you want to change the paramters, here is another example with the parameters we used for the chicken example in the manuscript.
   ````
   python clumppling -i input/chicken -o output/chicken_color -f structure --cd_param=1.05 --custom_cmap #D65859,#00AAC1,#01C0F6,#FDF0C4,#F1B38C,#AAD6BD,#6BB582,#B5DDF7,#AE8557,#FCEC73,#A4A569,#4264AC,#A1CDB2,#DE9D5D,#D9439A,#ABB2BA,#8775B3,#B3865C,#DADDE6,#E7BDD1,#FF9999"
   ````
   Here we are setting the community detection parameter ``--cd_param`` to be 1.05 and providing a custom colormap for visualization of the alignment results. The outputs will be saved in the directory ``output/chicken_color`` and a zipped file with the same contents and a same name will also be created.

* You may also put these command in a bash script ``ex2_brief.sh`` and run the script directly
   ````
   bash scripts/ex2_brief.sh
   ````
   More example commands to call the program can be found under the ```scripts``` directory.


**Detailed instructions to run the program can be found in the Manual.**

## Reference
Manuscript in prep.

