# ArboPhyl version 1.1
# 04-04-2023
# Author: Tim Verschuren

import subprocess
import os
import collections
from Bio import AlignIO
import shutil

def title():
    print(f"+-------------------------------------------------+")
    print(f"|\033[92m    ___       _\033[00m           ______ _           _   |")
    print(f"|\033[92m   / _ \     | |\033[00m          | ___ \ |         | |  |")
    print(f"|\033[92m  / /_\ \_ __| |__   ___ \033[00m | |_/ / |__  _   _| |  |")
    print(f"|\033[92m  |  _  | '__| '_ \ / _ \\\033[00m |  __/| '_ \| | | | |  |")
    print(f"|\033[92m  | | | | |  | |_) | (_) |\033[00m| |   | | | | |_| | |  |")
    print(f"|\033[92m  \_| |_/_|  |_.__/ \___/\033[00m \_|   |_| |_|\__, |_|  |")
    print(f"|                     ----------+        __/ |    |")
    print(f"|     --------------------+     |-------|___/     |")
    print(f"|                         |-----+                 |")
    print(f"|           --------------+                       |")
    print(f"+-------------------------------------------------+")
    print("\n")
    print("Welcome to \033[92mArbo\033[00mPhyl!\n| 0: Full pipeline\n| 1: BUSCO\n| 2: FilterBUSCOs\n| 3: MAFFT\n| 4: TrimAl\n| 5: MSA QC\n| 6: IQTREE Model Prediction\n| 7: Partition file creation\n| 8: IQTREE\n| q: Quit ArboPhyl\n")


def analysis_selection():
    while True:
        input_list = input("Select number(s) of required analyses seperated by commas (1,2,3): ")
        analyses = input_list.split(",")
        # Make sure that a program is not run multiple times
        if len(analyses) >> 1 & analyses.__contains__("0"):
            analyses = ["0"]
        # Allow user to exit program
        if analyses.__contains__("q"):
            quit()
        try:
            for item in analyses:
                int(item)
        except:
            print("\033[1;31mPlease select a valid input!\033[00m\n")
        else:
            break

    lineage = input("Please input the lineage used during the BUSCO analysis: ").lower()

    # Input number of threads, and check for correct input
    while True:
        try:
            print(f"\033[1;34mThere are {os.cpu_count()} threads available\033[00m")
            threads = input("Please input the number of threads to be used: ")
            int(threads)
        except KeyboardInterrupt:
            break
        except:
            print("\033[1;31mPlease only input integers!\033[00m\n")
        else:
            break

    return(analyses, lineage, threads)


def progress_bar(progress, total):
    percent = 100 * (progress / float(total))
    bar = '█' * int(percent) + '-' * (100 - int(percent))
    print(f"\r Progress: |{bar}| {percent:.0f}%", end = "\r")


def find_overlap(lineage):

    total_folders = 0
    busco_dict = {}
    folders = os.listdir()

    for file in folders:
        if file.__contains__("."):
            pass
        elif file.__contains__("busco_downloads"):
            pass
        elif file.__contains__("FilterBUSCOs_output"):
            pass
        else:
            total_folders += 1
            busco_list = []
            # Creates the filepath based on submitted database
            path = f"{os.path.dirname(os.path.abspath(file))}/{file}/run_{lineage}_odb10/busco_sequences/single_copy_busco_sequences"
            for file_name in os.listdir(path):
                # Place all fasta files into a list
                if file_name.endswith('.fna'):
                    busco_list.append(file_name)
                    if file_name not in busco_dict.keys():
                        busco_dict[file_name] = 1
                    else:
                        busco_dict[file_name] += 1
                else:
                    pass
    return [busco_dict, total_folders, lineage]

        
def busco_MSA(inputs):

    folders= os.listdir()
    i = 0

    for key, value in inputs[0].items():
        if value == inputs[1]:
            output_file_name = f"FilterBUSCOs_output/{key.replace('.fna', '_MS.fna')}"
            os.makedirs(os.path.dirname(output_file_name), exist_ok=True)
            with open(output_file_name, "w") as MSA_fasta:
                for file in folders:
                    if file.__contains__("."):
                        pass
                    elif file.__contains__("busco_downloads"):
                        pass
                    elif file.__contains__("FilterBUSCOs_output"):
                        pass
                    else:
                        path = f"{os.path.dirname(os.path.abspath(file))}/{file}/run_{inputs[2]}_odb10/busco_sequences/single_copy_busco_sequences"
                        for file_name in os.listdir(path):
                            if file_name == key:
                                with open(f"{path}/{file_name}", "r") as fasta_read:
                                    fasta_seq = fasta_read.readlines()
                                    for line in fasta_seq:
                                        if line.__contains__(">"):
                                            MSA_fasta.write(line.replace(line, f">{file}\n"))
                                        else:
                                            MSA_fasta.write(line)
            i += 1
            progress_bar(i, collections.Counter(inputs[0].values())[inputs[1]])
    print("\n")


def create_partition_file():

    model_dict = {}
    path = f"{os.getcwd()}/FilterBUSCOs_output/MAFFT_output/Trimmed_MSAs/Passed_MSA/Models"
    with open("Partition.nex", "w") as partition_file:
        partition_file.write("#nexus\nbegin sets;\n")
        for dir in os.listdir(path):
            if dir.__contains__(".") or dir.__contains__("_"):
                pass
            else:
                partition_file.write(f"\tcharset {dir} = FilterBUSCOs_output/MAFFT_output/Trimmed_MSAs/Passed_MSA/{dir}_trimmed.fna: *;\n")
                with open(f"{path}/{dir}/{dir}_trimmed.fna.iqtree", "r") as models:
                    for line in models:
                        if "BIC:" in line:
                            model_dict[dir] = (line.split(": ")[1]).replace("\n", "")
        partition_file.write("\tcharpartition mine = ")
        for key, value in model_dict.items():
            partition_file.write(f"{value}:{key}, ")
        partition_file.write(";\nend;")


def sum_of_pairs():

    path = f"{os.getcwd()}/FilterBUSCOs_output/MAFFT_output/Trimmed_MSAs"
    files = os.listdir(path)
    os.mkdir(f"{path}/Passed_MSA", )
    os.mkdir(f"{path}/Failed_MSA")

    for file in files:
        if file.endswith(".fna"):
            alignment = AlignIO.read(f"{path}/{file}", "fasta")
            sum_pair_score = 0

            for i in range(len(alignment[0])):
                pair_score = 0
                pair_list = []
                for j in range(len(alignment)):
                    pair_list.append(alignment[j].seq[i])
                score = 0
                for x in range(len(pair_list)):
                    for y in range(x+1, len(pair_list)):
                        if pair_list[x] == pair_list[y]:
                            score += 1
                pair_score = score/(len(pair_list)*(len(pair_list)-1)/2)
                sum_pair_score += pair_score
            if sum_pair_score/len(alignment[0]) >= 0.45:
                shutil.move(f"{path}/{file}", f"{path}/Passed_MSA/{file}")
            else:
                shutil.move(f"{path}/{file}", f"{path}/Failed_MSA/{file}")


def analysis_exe(inputs):
    # Start BUSCO
    if inputs[0].__contains__("0") or inputs[0].__contains__("1"):
        print(subprocess.run(["bash", "-i", "Arbophyl.sh", "busco", inputs[1], inputs[2]]))
    # Start FilterBUSCOs
    if inputs[0].__contains__("0") or inputs[0].__contains__("2"):
        busco_MSA((find_overlap(inputs[1])))
    # Start MAFFT
    if inputs[0].__contains__("0") or inputs[0].__contains__("3"):
        print(subprocess.run(["bash", "-i", "Arbophyl.sh", "mafft", inputs[2]]))
    # Start TrimAl
    if inputs[0].__contains__("0") or inputs[0].__contains__("4"):
        print(subprocess.run(["bash", "-i", "Arbophyl.sh", "trimal"]))
    # Start MSA QC
    if inputs[0].__contains__("0") or inputs[0].__contains__("5"):
        sum_of_pairs()
    # Start IQTREE Models
    if inputs[0].__contains__("0") or inputs[0].__contains__("6"):
        print(subprocess.run(["bash", "-i", "Arbophyl.sh", "iqtree_models", inputs[2]]))
    # Start Partition file creation
    if inputs[0].__contains__("0") or inputs[0].__contains__("7"):
        create_partition_file()
    # Start IQTREE
    if inputs[0].__contains__("0") or inputs[0].__contains__("8"):
        print(subprocess.run(["bash", "-i", "Arbophyl.sh", "iqtree", inputs[2]]))


def analysis_exe_v2(inputs):
    # Start BUSCO
    if inputs[0].__contains__("0") or inputs[0].__contains__("1"):
        os.system("""
            for file in *fasta
            do
                busco -i $file -o ${file/_*/""}/ -l $2 -m genome -c $3
            done
        """)
    # Start FilterBUSCOs
    if inputs[0].__contains__("0") or inputs[0].__contains__("2"):
        busco_MSA((find_overlap(inputs[1])))
    # Start MAFFT
    if inputs[0].__contains__("0") or inputs[0].__contains__("3"):
        print(subprocess.run(["bash", "-i", "Arbophyl.sh", "mafft", inputs[2]]))
    # Start TrimAl
    if inputs[0].__contains__("0") or inputs[0].__contains__("4"):
        print(subprocess.run(["bash", "-i", "Arbophyl.sh", "trimal"]))
    # Start MSA QC
    if inputs[0].__contains__("0") or inputs[0].__contains__("5"):
        sum_of_pairs()
    # Start IQTREE Models
    if inputs[0].__contains__("0") or inputs[0].__contains__("6"):
        print(subprocess.run(["bash", "-i", "Arbophyl.sh", "iqtree_models", inputs[2]]))
    # Start Partition file creation
    if inputs[0].__contains__("0") or inputs[0].__contains__("7"):
        create_partition_file()
    # Start IQTREE
    if inputs[0].__contains__("0") or inputs[0].__contains__("8"):
        print(subprocess.run(["bash", "-i", "Arbophyl.sh", "iqtree", inputs[2]]))


os.system("clear")
title()
analysis_exe(analysis_selection())

print("Finished all processes, thank you for using \033[92mArbo\033[00mPhyl!")
