#!/bin/bash

# Function generating progress bar of analyses
function progress_bar(){
	progress=$1
	total=$2
	declare -i percentage=100*$progress/$total
	declare -i _percentage=100-$percentage

	_done=$(printf "%${percentage}s")
	_todo=$(printf "%${_percentage}s")
	printf "\r Progress: |${_done// /█}${_todo// /-}| $percentage%%"
}

# Commands running BUSCO
if [[ "$1" == "busco" ]]; then
    mkdir $3
    mkdir $3/BUSCO_output/
    conda activate busco
    for file in $2*.f*
    do
        filename=${file##*/}
        cp $file $3/$filename
    done
    cd $3
    for file in *.f*
    do
        filename=${file##*/}
        busco -i $filename -o BUSCO_output/${filename/.*/""}/ -m $4 -l $5 -c $6
        rm $filename
    done
fi

# Commands running MAFFT
if [[ "$1" == "mafft" ]]; then
    cd $2
    conda activate mafft
    total_folder=$(ls -1 Filtered_BUSCOs/| wc -l)
    echo "Running MAFFT..."
    echo
    for file in Filtered_BUSCOs/*
    do
        filename=${file##*/}
        mafft --auto --inputorder --thread $3 --quiet "$file" > "MAFFT_output/${filename/MS/"MSA"}"
        next_folder=$(ls -1 MAFFT_output/ | wc -l)
        progress_bar $next_folder $total_folder
    done
    echo
    echo
fi


# Commands running trimAl
if [[ "$1" == "trimal" ]]; then
    mkdir $2/Models/
    cd $2
    conda activate trimal
    total_folder=$(ls -1 MAFFT_output/| wc -l)
    echo "Running trimAl..."
    echo
    for file in MAFFT_output/*
    do
        filename=${file##*/}
        trimal -in $file -out Models/${filename/MSA/"trimmed"} -strict > /dev/null 2>&1
        next_folder=$(ls -1 Models/ | wc -l)
        progress_bar $next_folder $total_folder
    done
    echo
    echo
fi

# Commands running IQTREE Model Prediction
if [[ "$1" == "iqtree_models" ]]; then
    cd $2
    for file in Models/*
    do
        filename=${file##*/}
        mkdir Models/${filename/_trimmed.f*a/""}
        mv $file Models/${filename/_trimmed.f*a/""}/$filename
    done

    cd Models/
    conda activate iqtree
    total_folder=$(ls -1 | wc -l)
    next_folder=0
    echo "Running IQTREE Model Prediction..."
    echo
    for dir in *
    do
        next_folder=$((next_folder + 1))
        cd $dir
        iqtree -s *.f* -m MF -nt $3 -quiet
        progress_bar $next_folder $total_folder
        cd ..
    done
    echo
    echo
fi

# Commands running IQTREE
if [[ "$1" == "iqtree" ]]; then
    cd $2
    conda activate iqtree
    iqtree -p *.nex -bb 1000 -alrt 1000 -nt $3
fi
