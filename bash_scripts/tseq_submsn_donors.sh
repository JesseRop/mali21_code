#!/bin/bash

##Variable order
#Pf chromosome, snp position, gene name, plasmoDB id, number of bases that flank snp, bam files input directory, plots output directory"
module load ISG/R/4.1.0

## Job submission variables
mem='-R "span[hosts=1] select[mem>160000] rusage[mem=160000]" -M 160000'
ncores=14

##job variables
list_item=($(seq 1 1 28))
# list_item=($(seq 3 1 8))
# list_item=(3 4 5 6 7 8 9 10 11 12 13 14)
# list_item=(1 2)
# list_item=(8)

exec_script=/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/tradeseq_donors.R

w_dir=/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/pseudotime_tradeseq_op

for f in "${!list_item[@]}";
do
    mkdir -p $w_dir/logs
	eval $(echo "bsub -q normal -n $ncores $mem -eo $w_dir/logs/tseq_job${list_item[f]}.err -oo $w_dir/logs/tseq_job${list_item[f]}.out Rscript $exec_script ${list_item[f]}")

done

