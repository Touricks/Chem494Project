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
3. Execute the program with command line:
```
python src/main.py
```
## Project architecture

```
ProjectFolder
└───data
    ├───prediction
        |————aptr6
        |————aptr9
        |————.....
    |————evaluation
        |————.....
    |————primer.txt
└───dataset
└───ctfiledatabase
└───src
    ├───Aptamer.py
    ├───Sequence.py
    ├───Kmer.py
    ├───Cluster.py
    ├───Round.py
    ├───main.py
```

## Parameter meaning
Parameter = [lenofseq, goal, num1, num2, op]

The primer.txt file in the data directory determines the starting position of N30. Note this!

lenofseq: The length of some sequence data is not 30 base pairs, and it can be adjusted here.

goal: The number of alternative aptamers. When the evaluation method is not used, the final recommended number of aptamers will be as per this value.

num1: The number of characteristic base short chains (kmers) counted when kmerlen < 8, recommended to be 50.

num2: The number of characteristic base short chains (kmers) counted when kmerlen > 8, recommended to be 1000.

    Note: The values of num1 and num2 will affect the efficiency of the algorithm.
    
op: Whether to allow the reuse of ct files. When op=True, the reuse of ct files is allowed, which can speed up the algorithm's efficiency.

Warning: The ct file contains basic information about the composition of the sequence. If changes are to be made to the sequence file itself (such as modifying lenofseq), op should be set to False, as the ct file is a file containing information about the hairpin structure of the sequence, which needs to be updated in this case.

## Prediction function intro
If you want to predict the aptamer for a substrate molecule based on a round of data, you just need to put the prefix of that round's data file name into pred=[].

## Evaluation function intro
Predictions based on a single round of data are inevitably biased, but we can use some basic methods to eliminate obviously incorrect predicted sequences.

Assume we make predictions based on the 9th round of apt data, which corresponds to setting pred=['aptr9'] in the program.

We can use the 5th round of apt data for preliminary screening, at which point eval=['aptr6', 'aptr9'] should be set.

If a predicted sequence has a higher abundance in the 6th round than in the 9th round, it is clearly not an aptamer. Such sequences can be removed in advance.

The remaining sequences will be saved in dataset/Recommendation after evaluation for target xxx.csv.

If the evaluation feature is not used, set eval=[]. In this case, the recommended sequences will be saved in dataset/Recommendation for target xxx.csv.
