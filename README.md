### Chem494Project
## Installation
This program is designed to run on a Linux environment since RNAFold repo is easy to install on Linux system.

Install all the packages in the requirements.txt file (all are available via pip).

1.1 Create virtual environment (venv) in the root directory of the project as shown in this [guide](https://docs.python.org/3/library/venv.html)

1.2 Execute the following commands in the terminal from the root directory of the project:
```
source venv/bin/activate
pip install -r requirements.txt
deactivate
```
2. Install RNAFold as indicated in this [guide](https://algosb2019.sciencesconf.org/data/RNAtutorial.pdf)

3. Check that the command RNAfold --help works in the terminal. If not, run the following commands before execution:

```
export PATH=${HOME}/Tutorial/Progs:${PATH}
```

## Execution instructions
1. Add data to data/ folder.
2. Place one fastq file in it, indicating we would provide prediction based on data from this file.

    Filename can be anything but only one file should exist in this folder.
3. Enter the src folder and execute the program
```
python3 main.py 15
```
The second parameter (e.g. 15) indicates the number of recommended aptamer you want to see. Default is 15.
## Project architecture

```
OptimalAptamerFinder
└───data
    ├───rounddata.fastq
└───src
    ├───Aptamer.py
    ├───Sequence.py
    ├───Kmer.py
    ├───Cluster.py
    ├───Round.py
    ├───main.py
└───README.md
└───requirements.txt
```
