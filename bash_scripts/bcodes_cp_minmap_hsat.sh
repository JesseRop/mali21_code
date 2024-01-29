#!/bin/bash

##copy filtered barcodes from cellranger to a folder structure that will allow running of souporcell with both hisat2 and minimap mappers

## Job submission variables
source_irods=("5736STDY11771535" "5736STDY11771536" "5736STDY11771544" "5736STDY11771546" "5736STDY11771545")
sample_name=("msc1" "msc3" "msc13" "msc14" "msc1272")

algn_nm=("hsat" "minmap")


w_dir=/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data

for i in "${!sample_name[@]}";
do
	for j in "${!algn_nm[@]}";
    do
	
    cp $w_dir/raw/${sample_name[i]}/${source_irods[i]}_DB52e_star/filtered_feature_bc_matrix/barcodes.tsv.gz $w_dir/processed/DB52e_star/${sample_name[i]}/bcodes/${algn_nm[j]}/
    gunzip -f $w_dir/processed/DB52e_star/${sample_name[i]}/bcodes/${algn_nm[j]}/barcodes.tsv.gz

    done
done


