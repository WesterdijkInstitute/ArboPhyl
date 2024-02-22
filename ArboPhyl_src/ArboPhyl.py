#!/usr/bin/env python

"""ArboPhyl is a pipeline for the automated construction of phylogenetic trees
based on BUSCO genes. 

usage: arbophyl [-h] -i INPUT -o OUTPUT -p PIPELINE 
[-m {genome,proteins,transcriptome}] [-l LINEAGE] [-s SHARED] [-c CPUS]
"""

"""Import Statements"""
import os
import re
import glob
from Bio import SeqIO


"""Authorship Information"""
__author__ = "Tim Verschuren"
__credits__ = ["Tim Verschuren", "Jérôme Collemare"]

__licence__ = "MIT"
__date__ = "22-02-2024"
__version__ = "1.0.0"
__maintainer__ = "Tim Verschuren"
__email__ = "t.verschuren@wi.knaw.nl"
__status__ = "Development"


def progress_bar(progress: int, total: int):
    """Progress bar displaying progress of ongoing analysis.

    Args:
        progress (int): Current increment.
        total (int): Total number of increments.
    """
    percent = 100 * (progress / float(total))
    bar = '█' * int(percent) + '-' * (100 - int(percent))
    print(f"\r Progress: |{bar}| {round(percent)}%", end = "\r")


def read_fasta(fasta_file) -> str:
    """Read content of fasta file and store sequence as string.

    Attributes:
        fasta_file (str): Path to the fasta file.

    Returns:
        seq (str): String containing sequence present in fasta file.
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()

    return seq


def auto_linebreak(string: str) -> str:
    """Automatically line break string at 60 characters.

    Args:
        string (str): String to be broken

    Returns:
        str: Line broken string
    """
    return re.sub("(.{60})", "\\1\n", string, 0, re.DOTALL)


class ap_analyses:
    """Class containing non-conda based analyses for ArboPhyl.
    """
    def __init__(self, input= "", output="", mode="", 
                 pipeline="", lineage="", shared=100) -> None:
        """Read the input settings of non-conda based analyses for ArboPhyl.

        Args:
            input (str, optional): File path to input folder.
            output (str, optional): File path to output folder.
            mode (str, optional): Analysis mode (genome, proteins).
            pipeline (list, optional): List containing to be executed 
            ArboPhyl modules
            lineage (str, optional): Name of lineage used in BUSCO analysis.
            shared (int, optional): Percentage of BUSCO genes that need to be 
            shared.
        """
        self.input = input
        self.output = output
        self.mode = mode
        self.pipeline = pipeline
        self.lineage = lineage
        self.shared = shared

    def get_BUSCOs(self) -> dict:
        """Retrieve paths of present BUSCO genes for each analysed species

        Returns:
            dict: Dictionary with species as the key and a list of paths 
            as its value
        """
        buscos_dict = {}
        print("Retrieving BUSCOs...\n")
        progress = 0

        #Look for nucleotide or amino acid sequences
        if self.mode == "genome":
            extention = ".fna"
        elif self.mode == "proteins":
            extention = ".faa"

        # Get path of BUSCO outputs and create list of folders
        path = f"{self.output}BUSCO_output/"
        folders = os.listdir(path)
        # Get list of present single copy busco genes per analysed species
        for folder in folders:
            if not folder.__contains__("Filtered_BUSCOs") and\
            not folder.__contains__("busco_downloads"):
                buscos_dict[folder] = []
                for file in os.listdir(glob.glob(
                    f"{path}/{folder}/*_odb10/busco*/single*")[0]):
                    if file.endswith(extention):
                        buscos_dict[folder].append(glob.glob(
                        f"{path}/{folder}/*_odb10/busco*/single*/{file}")[0])
            progress += 1
            progress_bar(progress, len(folders))
        print("\n")
        
        return buscos_dict

    @staticmethod
    def get_Overlap(busco_dict: dict) -> dict:
        """Calculate the percentage of species which possess specific 
        BUSCO genes.

        Args:
            busco_dict (dict): Output of get_BUSCOs, containing 
            species and paths to BUSCO genes

        Returns:
            dict: Dictionary containing the percentage of species which posses 
            a BUSCO gene as values and the BUSCO genes as keys.
        """
        overlap = {}
        # Count occurence of each BUSCO gene
        for value in busco_dict.values():
            for gene in value:
                if gene.split("/")[-1] in overlap:
                    overlap[gene.split("/")[-1]] += 1
                else:
                    overlap[gene.split("/")[-1]] = 1
        # Calculate percentage of species that possess BUSCO gene
        overlap = {k: round(v / len(busco_dict)*100, 1) \
                   for k, v in overlap.items()}
        return overlap

    def filter_BUSCOs(self, overlap: dict, buscos: dict):
        """Write overlapping BUSCOs to fasta files containing each species' 
        copy of one specific BUSCO gene

        Args:
            overlap (dict): Output of get_Overlap: Dictionary containing the 
            percentage of species which posses a BUSCO gene as values and the 
            BUSCO genes as keys. 
            buscos (dict): Output of get_BUSCOs: Dictionary with species as the 
            key and a list of paths as its value
        """
        progress = 0
        print("Filtering BUSCOs...\n")
        # Path to output locations
        path = f"{self.output}/Filtered_BUSCOs/"
        # Create output folder
        os.makedirs(os.path.dirname(path), exist_ok=True)

        # Per passed gene, create a fasta file containing the shared
        # sequences between each species
        for gene, value in overlap.items():
            present_dict = {}
            if value >= self.shared:
                # Determine analysis mode
                if self.mode == "genome":
                    MS_file = gene.split("/")[-1].replace('.fna', '_MS.fna')
                elif self.mode == "proteins":
                    MS_file = gene.split("/")[-1].replace('.faa', '_MS.faa')
                # Determine which species have fasta files available
                for species, busco_list in buscos.items():
                    if any(paths.endswith(gene) for paths in busco_list) \
                        == True:
                        present_dict[species] = True
                    else:
                        present_dict[species] = False

            # For each species, write shared gene sequence to file
            if len(present_dict) > 0:
                with open(f"{path}{MS_file}", "w") as MS_fasta:
                    for species, bool_val in present_dict.items():
                        if bool_val == True:
                            for paths in buscos[species]:
                                if paths.split("/")[-1] == gene:
                                    # Write sequence to file. Break line at 
                                    # 60 characters
                                    MS_fasta.write(f">{species}\n")
                                    MS_fasta.write(auto_linebreak(
                                                   read_fasta(paths)))
                                    MS_fasta.write("\n")
                        # If gene sequence not avaiable, replace with gaps
                        else:
                            for paths in buscos[list(present_dict.keys())\
                                    [list(present_dict.values()).index(True)]]:
                                if paths.split("/")[-1] == gene:
                                    gap = "-"*len(read_fasta(paths))
                                    # Write sequence to file. Break line at 
                                    # 60 characters
                                    MS_fasta.write(f">{species}\n")
                                    MS_fasta.write(auto_linebreak(gap))
                                    MS_fasta.write("\n")
            progress += 1
            progress_bar(progress, len(overlap))
        print("\n")

    def create_partition(self) -> None:
        """Read models from iqtree files and write to partition file.
        """
        progress = 0
        print("Creating partition file...\n")

        model_dict = {}
        path = f"{self.output}/Models/"

        # Generate name of partition file
        if self.output.endswith("/"):
            output_name = self.output.split("/")[-2]
        else:
            output_name = self.output.split("/")[-1]

        # Open partition file and write header line
        with open(f"{self.output}/{output_name}.nex", "w") as partition_file:
            partition_file.write("#nexus\nbegin sets;\n")
            # For each folder within Models, get names of files
            for dir in os.listdir(path):
                dir_name = dir.split("/")[-1]
                if self.mode == "genome":
                    file_name = f"{dir.split('/')[-1]}_trimmed.fna"
                elif self.mode == "proteins":
                    file_name = f"{dir.split('/')[-1]}_trimmed.faa"
                # Write path to file in partition file
                partition_file.write(
                    f"\tcharset {dir_name} = "\
                    f"Models/{file_name.split('_')[0]}/{file_name}: *;\n")
                # Read best model from .iqtree file and store in dictionary
                with open(f"{path}/{dir_name}/{file_name}.iqtree", 
                          "r") as models:
                    for line in models:
                        if "BIC:" in line:
                            model_dict[dir] = (line.split(": ")[1])\
                                .replace("\n", "")
                # Update progress bar
                progress += 1
                progress_bar(progress, len(os.listdir(path))*2)

            partition_file.write("\tcharpartition mine = ")
            # For each dictionary entry, write used model to partition file
            for key, value in model_dict.items():
                partition_file.write(f"{value}:{key}, ")
                progress += 1
                progress_bar(progress, len(os.listdir(path))*2)
            partition_file.write(";\nend;")
            print("\n")

