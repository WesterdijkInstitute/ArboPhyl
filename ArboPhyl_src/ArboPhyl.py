#!/usr/bin/env python

"""ArboPhyl is a pipeline for the automated construction of phylogenetic trees
based on BUSCO genes. 

usage: arbophyl [-h] -i INPUT -o OUTPUT -p PIPELINE 
[-m {genome,proteins,transcriptome}] [-l LINEAGE] [-s SHARED] [-c COMPLETE] 
[-t THREADS]
"""

"""Import Statements"""
import os
import re
import sys
import glob
from Bio import SeqIO


"""Authorship Information"""
__author__ = "Tim Verschuren"
__credits__ = ["Tim Verschuren", "Jérôme Collemare"]

__licence__ = "MIT"
__date__ = "07-03-2024"
__version__ = "1.1.0"
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


def read_fasta(fasta_file: str) -> str:
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
    def __init__(self, input= "", output="", mode="", pipeline="", 
                 lineage="", shared=100., complete=0.) -> None:
        """Read the input settings of non-conda based analyses for ArboPhyl.

        Args:
            input (str, optional): File path to input folder.
            output (str, optional): File path to output folder.
            mode (str, optional): Analysis mode (genome, proteins).
            pipeline (list, optional): List containing to be executed 
            ArboPhyl modules
            lineage (str, optional): Name of lineage used in BUSCO analysis.
            shared (float, optional): Percentage of BUSCO genes that need to be 
            shared.
            complete (float, optional): Minimum BUSCO completeness of genomes.
        """
        self.input = input
        self.output = output
        self.mode = mode
        self.pipeline = pipeline
        self.lineage = lineage
        self.shared = shared
        self.complete = complete

    def get_BUSCOs(self) -> dict:
        """Retrieve paths of present BUSCO genes for each analysed species

        Returns:
            dict: Dictionary with species as the key and a list of paths 
            as its value
        """
        buscos_dict = {}
        print("Retrieving BUSCOs...\n")

        #Look for nucleotide or amino acid sequences
        if self.mode == "genome":
            extention = ".fna"
        elif self.mode == "proteins":
            extention = ".faa"

        # Get path of BUSCO outputs and create list of folders
        path = f"{self.output}BUSCO_output/"
        folders = os.listdir(path)
        # Get list of present single copy busco genes per analysed species
        for progress, folder in enumerate(folders, 1):
            buscos_dict[folder] = []
            # For each file that matches wildcards, add filepath to dict
            for file in os.listdir(glob.glob(
                f"{path}/{folder}/*_odb10/busco*/single*")[0]):
                if file.endswith(extention):
                    buscos_dict[folder].append(glob.glob(
                    f"{path}/{folder}/*_odb10/busco*/single*/{file}")[0])
            progress_bar(progress, len(folders))
        print("\n")

        # Remove genomes with low completeness
        for key, val in self.BUSCO_qc().items():
            if val < self.complete:
                del buscos_dict[key]

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

    def filter_BUSCOs(self, overlap: dict, buscos: dict) -> None:
        """Write overlapping BUSCOs to fasta files containing each species' 
        copy of one specific BUSCO gene

        Args:
            overlap (dict): Output of get_Overlap: Dictionary containing the 
            percentage of species which posses a BUSCO gene as values and the 
            BUSCO genes as keys. 
            buscos (dict): Output of get_BUSCOs: Dictionary with species as the 
            key and a list of paths as its value
        """
        # Check if any gene matches the shared percentage, else exit analysis.
        if any(perc >= self.shared for perc in overlap.values()) == False:
            sys.exit(f"\033[1;31mWARNING: No BUSCOs detected that matched the"\
                        f" submitted shared percentage ({self.shared}%), "\
                        f"please lower shared percentage to at least "\
                        f"{max(overlap.values())}%.\033[00m\n")
        
        print("Filtering BUSCOs...\n")
        # Path to output locations
        path = f"{self.output}/Filtered_BUSCOs/"
        # Create output folder
        os.makedirs(os.path.dirname(path), exist_ok=True)

        # Per passed gene, create a fasta file containing the shared
        # sequences between each species
        for progress, (gene, value) in enumerate(overlap.items(), 1):
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
            progress_bar(progress, len(overlap))
        print("\n")

    def create_partition(self) -> None:
        """Read models from iqtree files and write to partition file.
        """
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
            for progress, dir in enumerate(os.listdir(path), 1):
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
                progress_bar(progress, len(os.listdir(path))*2)

            partition_file.write("\tcharpartition mine = ")

            # For each dictionary entry, write used model to partition file
            for progress, (key, value) in enumerate(model_dict.items(), 
                                                    progress+1):
                partition_file.write(f"{value}:{key}, ")
                progress_bar(progress, len(os.listdir(path))*2)
            partition_file.write(";\nend;")
            print("\n")

    def BUSCO_qc(self) -> dict:
        """_summary_

        Returns:
            dict: _description_
        """
        comp_dict = {}

        # Get path of BUSCO outputs and create list of folders
        path = f"{self.output}BUSCO_output/"
        folders = os.listdir(path)
        # Retrieve BUSCO completeness from analyses
        for folder in folders:
            for file in glob.glob(f"{path}/{folder}/*txt"):
                with open(file, "r") as summary:
                    completeness = summary.readlines()[8].split(":")[1]\
                        .split("%")[0]
                    comp_dict[folder] = float(completeness)

        return comp_dict

    def BUSCO_qc_screen(self, qc_dict: dict) -> None:
        """Prints a window showing BUSCO completeness of genomes.

        Args:
            qc_dict (dict): Dictionary containing completeness scores of each
            genome.
        """
        remove_list = []

        # Assign colour codes
        colours = {
            "g": ["\033[92m", "\033[00m"],
            "y": ["\033[93m", "\033[00m"],
            "r": ["\033[91m", "\033[00m"]
        }

        # Set thresholds to defaults or submitted value
        threshold = [95, 90] if self.complete == 0. \
                    else [self.complete, self.complete]

        # Determine the widths of the window frame
        width = len(max(qc_dict, key=len))
        frame_width = width + len(max([str(value) for value in \
                                       qc_dict.values()], key=len)) + 9
        
        # Print header
        print(f"{' '*((frame_width-28)//2)}BUSCO COMPLETENESS OF GENOMES")

        # Print top part of window
        print(f"+{'-'*(frame_width)}+")
        # For each genome, determine if the quality score reaches a threshold
        # and colour accordingly.
        for key, val in qc_dict.items():
            if val >= threshold[0]:
                print(f"|  {colours['g'][0]}{key}"\
                        f"{4*' ' + ' '*(width-len(key))}"\
                        f"{val}%{colours['g'][1]}  |")
            if val >= threshold[1] and val < threshold[0]:
                print(f"|  {colours['y'][0]}{key}"\
                        f"{4*' ' + ' '*(width-len(key))}"\
                        f"{val}%{colours['y'][1]}  |")
            if val < threshold[1]:
                if self.complete != 0.:
                    remove_list.append(key)
                print(f"|  {colours['r'][0]}{key}"\
                        f"{4*' ' + ' '*(width-len(key))}"\
                        f"{val}%{colours['r'][1]}  |")
        # Lower part of window
        print(f"+{'-'*(frame_width)}+\n")

        # If files need to be skipped, print names.
        if self.complete != 0. and len(remove_list) > 0:
            print(f"Skipping following genomes: {str(remove_list)[1:-1]}\n")
        elif self.complete == 0.:
            print(f"'Complete' option set to default, keeping all genomes.\n")

        