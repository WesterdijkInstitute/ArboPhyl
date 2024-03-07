#!/usr/bin/env python

"""Argument parsing and analysis initialization script for ArboPhyl. 
"""

"""Import Statements"""
from .ArboPhyl import *
import subprocess
import argparse
from argparse import RawTextHelpFormatter

def main():
    """Function executing ArboPhyl scripts.
    """
    argParser = argparse.ArgumentParser(description=
                            "ArboPhyl is a BUSCO based "\
                            "pipeline for the construction of " \
                            "phylogenetic trees.\n"\
                            "For more information see: "\
                            "https://github.com/WesterdijkInstitute/ArboPhyl",
                            formatter_class=RawTextHelpFormatter)
    argParser.add_argument("-i",
                           "--input",
                           type=str,
                           help="Path to input folder",
                           required=True)
    argParser.add_argument("-o",
                           "--output",
                           type=str,
                           help="Path to output location",
                           required=True)
    argParser.add_argument("-p",
                           "--pipeline",
                           type=str,
                           help="Processes to be performed in the pipeline, "\
                            "separated by commas (e.g., 1,2,3)\n"\
                            "0: Full pipeline\n"\
                            "1: BUSCO\n"\
                            "2: Filter BUSCOs\n"\
                            "3: MAFFT\n"\
                            "4: TrimAl\n"\
                            "5: IQTREE Model Prediction\n"\
                            "6: Partition file creation\n"\
                            "7: IQTREE",
                           required=True)
    argParser.add_argument("-m",
                           "--mode",
                           type=str,
                           help="Select mode based on desired input",
                           choices=["genome", "proteins"],
                           required=True)
    argParser.add_argument("-l",
                           "--lineage",
                           type=str,
                           help="BUSCO lineage (e.g., ascomycota)",
                           required=False)
    argParser.add_argument("-s",
                           "--shared",
                           type=float,
                           help="Percentage of BUSCOs that must be shared "\
                           "across analysed species (default: 100%%)",
                           required=False)
    argParser.add_argument("-c",
                           "--complete",
                           type=float,
                           help="Required BUSCO completeness of genomes. "\
                            "Keeps all sequences by default unless specified "\
                            "otherwise (e.g., 98%%)",
                           required=False)
    argParser.add_argument("-t",
                           "--threads",
                           type=str,
                           help="Number of CPUs for analyses, default: auto",
                           required=False)

    args = argParser.parse_args()
    # Default parameter for shared BUSCOs & BUSCO completeness
    if not args.shared:
        args.shared = 100.
    if not args.complete:
        args.complete = 0.
    
    # Make sure each analysis runs only once, and split input
    args.pipeline = args.pipeline.split(",")
    if args.pipeline.__contains__("0"):
        args.pipeline = "0"
    # Make sure no letters are submitted
    try:
        for item in args.pipeline:
            int(item)
    except:
        print("\033[1;31mPlease select a valid pipeline option!\033[00m\n")
        exit()

    # BUSCO parameters are required
    if not args.lineage:
        if  args.pipeline.__contains__("0") \
            or args.pipeline.__contains__("1") \
            or args.pipeline.__contains__("6"):
            print("\033[93m!!! BUSCO related parameter required " \
                  "(lineage) !!!\033[00m")
            exit()
    else:
        # Remove case sensitivity
        args.lineage = args.lineage.lower()

    # Shared percentage is required
    if not args.shared:
        if args.pipeline.__contains__("0") \
        or args.pipeline.__contains__("2"):
            print("\033[93m!!! Shared percentage required !!!\033[00m")
            exit()            

    # Run modules submitted in pipeline
    # Run BUSCO
    if args.pipeline.__contains__("0") or args.pipeline.__contains__("1"):
        if not args.threads:
            args.threads = "1"
        print(subprocess.run(["bash", "-i", "ArboPhyl.sh", "busco", 
                              args.input, args.output, args.mode, 
                              args.lineage, args.threads]))
    
    # Run Filter BUSCOs
    if args.pipeline.__contains__("0") or args.pipeline.__contains__("2"):
        busco_dict = ap_analyses(output=args.output, 
                                mode=args.mode,
                                complete=args.complete).get_BUSCOs()
        ap_analyses(complete=args.complete)\
            .BUSCO_qc_screen(ap_analyses(output=args.output).BUSCO_qc())
        ap_analyses(output=args.output, shared=args.shared, mode=args.mode)\
            .filter_BUSCOs(ap_analyses.get_Overlap(busco_dict), busco_dict)
        
    # Run MAFFT Multiple sequence alignment
    if args.pipeline.__contains__("0") or args.pipeline.__contains__("3"):
        if not args.threads:
            args.threads = "-1"
        os.makedirs(os.path.dirname(f"{args.output}/MAFFT_output/"), 
                    exist_ok=True)
        print(subprocess.run(["bash", "-i", "ArboPhyl.sh", "mafft", 
                              args.output, args.threads]))
    
    # Run TrimAl on MSAs
    if args.pipeline.__contains__("0") or args.pipeline.__contains__("4"):
        print(subprocess.run(["bash", "-i", "ArboPhyl.sh", "trimal", 
                              args.output]))

    # Run IQTREE model finder    
    if args.pipeline.__contains__("0") or args.pipeline.__contains__("5"):
        if not args.threads:
            args.threads = "AUTO"
        print(subprocess.run(["bash", "-i", "ArboPhyl.sh", "iqtree_models",
                              args.output, args.threads]))
    
    # Create partition file
    if args.pipeline.__contains__("0") or args.pipeline.__contains__("6"):
        ap_analyses(output=args.output, mode=args.mode).create_partition()

    # Run IQTREE 
    if args.pipeline.__contains__("0") or args.pipeline.__contains__("7"):
        if not args.threads:
            args.threads = "AUTO"
        print(subprocess.run(["bash", "-i", "ArboPhyl.sh", "iqtree",
                              args.output, args.threads]))
    
    title()
    print("Finished all processes, thank you for using "\
          "\033[92mArbo\033[00mPhyl!\n")


def title():
    """Title card of ArboPhyl.
    """
    print(f"+-------------------------------------------------+")
    print(f"|\033[92m    ___       _\033[00m           ______ _ "\
          "          _   |")
    print(f"|\033[92m   / _ \     | |\033[00m          | ___ \ |"\
          "         | |  |")
    print(f"|\033[92m  / /_\ \_ __| |__   ___ \033[00m | |_/ / |"\
          "__  _   _| |  |")
    print(f"|\033[92m  |  _  | '__| '_ \ / _ \\\033[00m |  __/| "\
          "'_ \| | | | |  |")
    print(f"|\033[92m  | | | | |  | |_) | (_) |\033[00m| |   | |"\
          " | | |_| | |  |")
    print(f"|\033[92m  \_| |_/_|  |_.__/ \___/\033[00m \_|   |_|"\
          " |_|\__, |_|  |")
    print(f"|                     ----------+        __/ |    |")
    print(f"|     --------------------+     |-------|___/     |")
    print(f"|                         |-----+                 |")
    print(f"|           --------------+                       |")
    print(f"+-------------------------------------------------+")
    print("\n")


if __name__ == "__main__":
    main()

