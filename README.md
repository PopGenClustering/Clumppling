# Clumppling

This is the GitHub repository for the program ***Clumppling*** (CLUster Matching and Permutation Program that uses integer Linear ProgrammING), a new framework for aligning clustering results of population structure analysis.

> The current version of *Clumppling* is **test version 0.0.1**.

## How to Install
1. Install [Python3](https://www.python.org/downloads/). Versions above 3.7 are recommended. If you use Git for installation, also make sure to install [Git](https://git-scm.com/).
2. Install the *Clumppling* package. This can be done by either one of the following two ways:
    * Download the package files at https://github.com/PopGenClustering/Clumppling/archive/refs/heads/master.zip and unzip to a local directory. Rename the folder named ``Clumppling-master`` to ``Clumppling``.
      
      Or, directly clone this repository to a local directory via ``git clone https://github.com/PopGenClustering/Clumppling.git``.
      
      Next, navigate into the directory ``Clumppling``, then run the following command
      ````
      pip install -e .
      ````

    * (Not recommended) Install the package directly from GitHub via ``pip install git+https://github.com/PopGenClustering/Clumppling``. Doing it this way will not get you all the example files automatically. To run examples, you still need to download the example files from this repository. 
  3. Check the installation by running
      ````
      python -m clumppling.main -h
      ````
      If the installation was successful, you should be prompted by the usage of this function.
    
## How to Run 
To run the program, follow instructions in the Manual.

## Example: Inferred Population Structure From Chicken Microsatellite Data

## Reference
Manuscript in prep.

