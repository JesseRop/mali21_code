#!/usr/bin/env nextflow

// 2022 cellranger output
// params.bam = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf/*/outs/possorted_genome_bam.bam" 

// 2021 cellranger output
params.bam = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/raw/msc*/*_DB52e_star/possorted_genome_bam.bam"

// params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/msc*/bcodes/highQ_cells.tsv" 
// params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/raw/msc*/*_DB52e_star/filtered_feature_bc_matrix/barcodes.tsv.gz"


// barcodes with doublets removed
// params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/msc14/bcodes/*/barcodes.tsv"
params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/msc*/bcodes/*/bcodes_no_dbt/barcodes.tsv"

params.o_dir= "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star"

ncores="20"
mem="140 GB"

params.soup_dir = "soupc_no_dbt"

params.scrpt = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/soupc_v25_hsat.sh"
params.ref  = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/hisat_refs/Pfalciparum.genome.fasta"

bam_ch = Channel
				.fromPath(params.bam)
				.map { file -> tuple(file.getParent().toString().split('\\/')[12], file, file+'.bai') }

// bam_ch.view()

bcodes_al_ch = Channel
				.fromPath(params.bcodes)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file.getParent().toString().split('\\/')[15], file) }
                .branch {
                    hsat_paths: it[1] == 'hsat'
                        return tuple(it[0], it[1], "HISAT2", it[2])
                    minmap_paths: it[1] == 'minmap'
                        return tuple(it[0], it[1], "minimap2", it[2])
                }

bcodes_ch = bcodes_al_ch.hsat_paths.concat(bcodes_al_ch.minmap_paths)

bcodes_bam_ch = bcodes_ch.combine(bam_ch, by:0)
// bcodes_bam_ch.view()

k_ch = Channel.from( 1..15 )
// k_ch.view()

process SOUPC1 {
    memory "${mem}"
    cpus "${ncores}"

    tag "Souporcell k1 ${sample_nm}  reads"

    // publishDir params.outdir_index_temp, mode: "copy"
    publishDir "$params.o_dir/${sample_nm}/${params.soup_dir}/${algn_nm}/parent/", mode: 'copy'

    input:
    tuple val(sample_nm), val(algn_nm), val(algnr), path(bcodes), path(bam), path(bai)
    // val(algnr)
	
    output:
    tuple val(sample_nm), val(algn_nm), path(out_dir)

	script:
	"""
    
    ${params.scrpt} ${bam} ${bcodes} 1 ${algnr} ${params.ref} out_dir ${ncores}
    
	"""			
    
}

// cp ${out_dir}/{vartrix.done,alt.mtx,ref.mtx,barcodes.tsv,bcftools.err,souporcell_merged_sorted_vcf.vcf.gz,souporcell_merged_sorted_vcf.vcf.gz.tbi,retagging.done,souporcell_minimap_tagged_sorted.bam.bai,souporcell_minimap_tagged_sorted.bam,retag.err,remapping.done,minimap.err,fastqs.done} out_dir_k${k}/
// cp ${bcodes} out_dir_k${k}/barcodes.tsv

process SOUPC_2PLUS {
    memory '60 GB'
    cpus '5'
    errorStrategy 'ignore'

    tag "Souporcell k${k} ${sample_nm}  reads"

    publishDir "$params.o_dir/${sample_nm}/${params.soup_dir}/${algn_nm}/k${k}/", mode: 'copy', overwrite: true

    input:
	tuple val(sample_nm), val(algn_nm), val(algnr), path(bcodes), path(bam), path(bai), path(out_dir)
    // val(algnr)
    each(k)
		
    output:
    // tuple val(sample_nm), path('*')
    // tuple val(sample_nm), path('out_dir_k*/clusters.tsv'), path('out_dir_k*/cluster_genotypes.vcf'), path('out_dir_k*/ambient_rna.txt'), path('clusters_log_likelihoods*')
    tuple val(sample_nm), val(algn_nm), path('clusters.tsv'), path('cluster_genotypes.vcf'), path('ambient_rna.txt'), path('clusters_log_likelihoods*')  

	script:
	"""
    export PATH="\$PATH:/software/team222/jr35/souporcell_hard_inst/souporcell/souporcell/target/release"
    export PATH="\$PATH:/software/team222/jr35/souporcell_hard_inst/souporcell/troublet/target/release"
    export PATH="\$PATH:/software/team222/jr35/souporcell_hard_inst/souporcell"


    ## Souporcell on vcf
    souporcell -a ${out_dir}/alt.mtx -r ${out_dir}/ref.mtx -b ${bcodes} -k ${k} -t ${task.cpus} > clusters_tmp.tsv 2>>clusters.err

    ## Doublet detction
    troublet -a ${out_dir}/alt.mtx -r ${out_dir}/ref.mtx --clusters clusters_tmp.tsv > clusters.tsv

    ## Genotype and ambient RNA
    consensus.py -c clusters.tsv -a ${out_dir}/alt.mtx -r ${out_dir}/ref.mtx --soup_out ambient_rna.txt -v ${out_dir}/souporcell_merged_sorted_vcf.vcf.gz -p 1 --vcf_out cluster_genotypes.vcf --output_dir .

    
    grep -H 'best total log probability' clusters.err > clusters_log_likelihoods_${k}.txt
	"""			
}


process LL_KNEE_PLOT {
    memory '10 GB'
    cpus '1'

    tag "Souporcell log likelihood values ${sample_nm} reads"

    publishDir "$params.o_dir/${sample_nm}/${params.soup_dir}/${algn_nm}/", mode: 'copy', overwrite: true

    input:
	tuple val(sample_nm), val(algn_nm), path(clst_err)
		
    output:
    // tuple val(sample_nm), path('*')
    tuple val(sample_nm), path('clusters_log_likelihoods.txt')
    
	script:
	"""
    cat $clst_err > clusters_log_likelihoods.txt
    
	"""			
}

workflow {
    
    soupc1_ch = SOUPC1(bcodes_bam_ch)
    soupc1_ch.view()

    bcodes_bam_soupc_ch = bcodes_bam_ch.combine(soupc1_ch, by:[0,1])
    bcodes_bam_soupc_ch.view()
    
    soupc2plus_ch = SOUPC_2PLUS(bcodes_bam_soupc_ch, k_ch)
    // soupc2plus_ch.view()
    soupc2plus_ll_ch = soupc2plus_ch
        // .map { file -> tuple(file[0], file[4]), file[4].getParent().getName() }
        .map { file -> tuple(file[0], file[1], file[5]) }
        .groupTuple(by: [0,1])
    
    soupc2plus_ll_ch.view()

    LL_KNEE_PLOT(soupc2plus_ll_ch)

}

