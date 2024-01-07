#!/usr/bin/env nextflow

// 2022 cellranger output
// params.bam = "/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/cellranger_runs/Pf/*/outs/possorted_genome_bam.bam" 

// 2021 cellranger output
// params.bam = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/raw/msc*/*_DB52e_star/possorted_genome_bam.bam"

// 2021 hisat2 bam
params.bam = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/msc14/soupc/minmap/parent/out_dir/souporcell_minimap_tagged_sorted.bam"
// params.bam = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/msc14/soupc/hsat/parent/out_dir/souporcell_minimap_tagged_sorted.bam"

// params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/msc*/bcodes/stage_afm_strain_k*/*.tsv" 
// params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/msc*_5736STDY*/bcodes/*strain_k*/*.tsv" 
// params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/msc*/bcodes/*strain_k*/*.tsv" 
// params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/msc*/bcodes/minmap/*strain_k*/*.tsv" 
// params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/msc14/bcodes/hsat/*strain_k{6,7}/*.tsv" 
params.bcodes = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/msc14/bcodes/minmap/*strain_qcd*/*.tsv"

params.odir= "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star"

// ncores=10
params.mapper="minmap"

params.scrpt = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/multipurpose_scripts/cr_subset_bam_linux.sh"
// params.ref  = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/Pfalciparum.genome.fasta"
params.ref  = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/hisat_refs/Pfalciparum.genome.fasta"

bam_ch = Channel
				.fromPath(params.bam)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file, file+'.bai') }

// bam_ch.view()

bcodes_ch = Channel
				.fromPath(params.bcodes)
				.map { file -> tuple(file.getParent().toString().split('\\/')[13], file.getParent().name, file.baseName, file) }

// bcodes_ch.view()
bcodes_bam_ch = bcodes_ch.combine(bam_ch, by:0)
bcodes_bam_ch.view()

process SUBSET {
    memory '150 GB'
    cpus '14'

    tag "Subset ${strn_stg}  reads"
    
    input:
    tuple val(sample_nm), val(grp), val(strn_stg), path(bcodes), path(bam), path(bai)
	path scrpt
	
    output:
    tuple val(sample_nm), val(grp), val(strn_stg), path(bcodes), path('*_sset_sorted.bam'), path('*_sset_sorted.bam.bai')

	script:
	"""
    ${params.scrpt} ${strn_stg} ${bcodes} ./ ${bam}
    
	"""			
    
}

// Process each group and ploidy
process TAGBAM {
    memory '100 GB'
    cpus '3'

    tag "Tagging ${strn_stg}  reads"

    // publishDir "$params.odir/${sample_nm}/pbulk_gtypes/${grp}_${params.mapper}/nw_bams", mode: 'copy', pattern: "*sset_sorted_rg.bam*"
        
    input:
    tuple val(sample_nm), val(grp), val(strn_stg), path(bcodes), path(bam), path(bai)

    output:
    tuple val(sample_nm), val(grp), val(strn_stg), path(bcodes), path('*_sset_sorted_rg.bam'), path('*_sset_sorted_rg.bam.bai')
    
    script:
    """
    module load common-apps/samtools/1.17

    samtools addreplacerg -r ID:${strn_stg} -r SM:${strn_stg} $bam -o ${strn_stg}_sset_sorted_rg.bam
    samtools index ${strn_stg}_sset_sorted_rg.bam
    """
}

// variant calling variables
// params.ref  = "/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/Pfalciparum.genome.fasta"

ref_ch = Channel.value(params.ref)
ploidy=[1]

// freebayes parameters
// fb_spec=Channel.of("-iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --min-alternate-fraction 0.2" "-C 2 -q 20 -n 3 -E 3 -m 30 --min-coverage 6 --min-alternate-fraction 0.2 --theta 0.01")
fb_spec=Channel.value("-iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --min-alternate-fraction 0.2")

// Bcftools vcf processing parameters
// bcf_spec=Channel.of("view --max-alleles 2" "norm --multiallelics -")
bcf_spec=Channel.value("view --max-alleles 2")


// Section 1##########################################################
// Run FreeBayes with the specified parameters

// Use ploidy 1 or 2 - souporcell uses 2 I think. We expect most cells to be homozygous - haploid
// ploidy=[1 2]

process FREEBAYES {
    memory '150 GB'
    cpus '28'

    tag "Freebayes variant calling on ${grp} ${sample_nm}"

    // publishDir "$params.outd", mode: 'copy'
    publishDir "$params.odir/${sample_nm}/pbulk_gtypes/${grp}_${params.mapper}/p${pld}/", mode: 'copy', pattern: "*.vcf", overwrite: true
    
    input:
    tuple val(sample_nm), val(grp), path(bam), path(bai)

    each pld
    val fb_s
    val bcf_s
    // val ref

    output:
    tuple val(sample_nm), val(grp), path('*_biallelic.vcf'), path('*_multiallelic.vcf')

    script:
    """
    module load common-apps/bcftools/1.17 
    module load freebayes/1.3.6--h346b5cb_1

    freebayes -f $params.ref $fb_s --ploidy $pld $bam > fb_multiallelic.vcf 
    bcftools $bcf_s fb_multiallelic.vcf -o fb_biallelic.vcf
    """
}

process VARTRIX_UMI {
    memory '100 GB'
    cpus '15'
    debug true

    tag "Vartrix allele counting on ${strn_stg} from ${grp} ${sample_nm}"
    
    // publishDir "$params.outd/vartrix_biallelic/", mode: 'copy'// 
    publishDir "$params.odir/${sample_nm}/pbulk_gtypes/${grp}_${params.mapper}/p${pld}/vartrix_biallelic", mode: 'copy', pattern: "*", overwrite: true
    
    input:
    tuple val(sample_nm), val(grp), val(strn_stg), path(bcodes), path(bam), path(bai), path(bi_vcf)
    each pld
    
    output:
    tuple val(sample_nm), val(grp), val(strn_stg), path ('*')
    // stdout
    
    script:
    """
    export PATH="\$PATH:/software/team222/jr35/vartrix/vartrix-1.1.22/target/release"
    vartrix --umi --out-variants ${strn_stg}_variants.txt --mapq 30 -b ${bam} -c ${bcodes} --scoring-method coverage --ref-matrix ${strn_stg}_ref.mtx --out-matrix ${strn_stg}_alt.mtx -v ${bi_vcf} --fasta ${params.ref}
    
    """
}


workflow {
    
    sset_bams_ch = SUBSET(bcodes_bam_ch, params.scrpt)
    // sset_bams_ch.view()
    tagd_bams_ch = TAGBAM(sset_bams_ch)
    // tagd_bams_ch.view()

    tagbam_ch_all = tagd_bams_ch
                        .map { file -> tuple(file[0], file[1], file[4], file[5]) }
                        .groupTuple(by: [0, 1])
    
    // tagbam_ch_all.view()

    fb_vcf_ch = FREEBAYES(tagbam_ch_all, ploidy, fb_spec, bcf_spec)

    // fb_vcf_ch.view()
    fb_vcf_bi_ch = fb_vcf_ch.map { file -> tuple(file[0], file[1], file[2]) }
    fb_vcf_bi_ch.view()

    bcodes_fb_vcf_ch = tagd_bams_ch
                        .combine(fb_vcf_bi_ch, by: [0, 1])

    bcodes_fb_vcf_ch.view()

    vartx_ch = VARTRIX_UMI(bcodes_fb_vcf_ch, ploidy)

    vartx_ch.view()

}
