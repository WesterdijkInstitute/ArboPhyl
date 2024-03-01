# ArboPhyl

ArboPhyl is a BUSCO based phylogenomics pipeline for the construction of phylogenetic trees from either nucleotide or protein fasta files. 

## Installation

***Python 3.10 is required for installation!***

ArboPhyl can be installed by first cloning this GitHub repository. Subsequently, go into the location where the repository was saved and use pip to install the package.
```bash
# Clone repository
git clone https://github.com/WesterdijkInstitute/ArboPhyl.git

# Go into repository location
cd .../ArboPhyl/

# Install package with pip
pip install .
```

The pipeline can be run on a local or server based linux environment with the [Anaconda](https://anaconda.org/) package manager installed. Once installed, the following tools need to be installed in their own environment:
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

ArboPhyl has four required inputs and three optional ones. The path to the input folder (containing either protein or nucleotide fasta files) as well as the output folder are required, in addition to the mode (genome/proteins) and the segments of the pipeline which need to be executed. The lineage input is only required when performing busco related analyses (e.g., 0, 1 and 2), the shared parameter - percentage of BUSCO genes that need to be shared across all analysed species - is set to 100% by default, but can be lowered. The cpus parameter uses the automatic settings for each analysis by default but can be specified by the user as well. 

Note: ***It is important that the input and output folders remain the same if the pipeline is executed in multiple segments instead of all at once.***
```
usage: arbophyl [-h] -i INPUT -o OUTPUT -p PIPELINE -m {genome,proteins} [-l LINEAGE] [-s SHARED] [-c CPUS]

ArboPhyl is a BUSCO based pipeline for the construction of phylogenetic trees.
For more information see: https://github.com/WesterdijkInstitute/ArboPhyl

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to input folder
  -o OUTPUT, --output OUTPUT
                        Path to output location
  -p PIPELINE, --pipeline PIPELINE
                        Processes to be performed in the pipeline, separated by commas (e.g., 1,2,3)
                        0: Full pipeline
                        1: BUSCO
                        2: Filter BUSCOs
                        3: MAFFT
                        4: TrimAl
                        5: IQTREE Model Prediction
                        6: Partition file creation
                        7: IQTREE
  -m {genome,proteins}, --mode {genome,proteins}
                        Select mode based on desired input
  -l LINEAGE, --lineage LINEAGE
                        BUSCO lineage (e.g., ascomycota)
  -s SHARED, --shared SHARED
                        Percentage of BUSCOs that must be shared across analysed species (default: 100%)
  -c CPUS, --cpus CPUS  Number of CPUs for analyses, default: auto
```

## Author
Tim Verschuren <br/>
[GitHub](https://github.com/TimVerschuren)
[LinkedIn](https://www.linkedin.com/in/tim-verschuren-27082919b/)
