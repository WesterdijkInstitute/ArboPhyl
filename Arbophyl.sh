#!/bin/bash

if [[ "$1" == "busco" ]]; then
    conda activate busco
    for file in *fasta
    do
        busco -i $file -o ${file/_*/""}/ -l $2 -m genome -c $3
    done
fi

if [[ "$1" == "mafft" ]]; then
    cd FilterBUSCOs_output/
    conda activate mafft
    mkdir MAFFT_output
    for file in *fna
    do
        mafft --auto --inputorder --thread $2 "$file" > "MAFFT_output/${file/MS/"MSA"}"
    done
fi

if [[ "$1" == "trimal" ]]; then
    cd FilterBUSCOs_output/MAFFT_output/
    conda activate trimal
    mkdir Trimmed_MSAs
    for file in *fna
    do
        trimal -in $file -out Trimmed_MSAs/${file/MSA/"trimmed"} -strict
    done
fi

if [[ "$1" == "iqtree_models" ]]; then
    cd FilterBUSCOs_output/MAFFT_output/Trimmed_MSAs/Passed_MSA
    mkdir Models
    for file in *fna
    do
        mkdir Models/${file/_trimmed.fna/""}
        cp $file Models/${file/_trimmed.fna/""}
    done

    cd Models/
    conda activate iqtree
    for dir in *
    do
    cd $dir
    iqtree -s *.fna -m MF -nt $2
    cd ..
    done
fi

if [[ "$1" == "iqtree" ]]; then
    conda activate iqtree
    iqtree -s FilterBUSCOs_output/MAFFT_output/Trimmed_MSAs/Passed_MSA/ -p Partition.nex -bb 1000 -alrt 1000 -nt $2 
fi