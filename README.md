# ArboPhyl

ArboPhyl is a pipeline for building phylogenetic trees for the placement of species into their correct position in the fungal tree of life. It uses overlapping BUSCO genes between various species to determine their differences and builds a tree using IQTREE. The goal of this pipeline is to facilitate the process of running the multiple analyses required to build a tree automatically with only a starting input by the user. 

## Installation

The pipeline can be run on a linux server with the [Anaconda](https://anaconda.org/) package manager installed. Once installed, the following tools need to be installed in their own environment:

```bash
# install BUSCO
conda create -n busco -c conda-forge -c bioconda busco=5.4.4

# install MAFFT
conda create -n mafft -c bioconda mafft

# install TrimAl
conda create -n trimal -c bioconda trimal

# install IQTREE
conda create -n iqtree -c bioconda iqtree
```

## Usage

Download the assembly scaffolds of each of your species and place their fasta files into one folder. Place both the [Arbophyl.py](https://github.com/TimVerschuren/ArboPhyl/blob/master/Arbophyl.py) and [Arbophyl.sh](https://github.com/TimVerschuren/ArboPhyl/blob/master/Arbophyl.sh) files into the folder aswell. Go to the correct directory utilising the commandline and run the following command:

```bash
python Arbophyl.py
```

Next select the analyses you want to carry out and the BUSCO database that is/was used by the BUSCO analysis.
