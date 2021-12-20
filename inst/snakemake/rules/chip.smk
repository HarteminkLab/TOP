
# Count ChIP-seq read coverage for candidate sites
# this requires getting the idxstats for (all) the ChIP-seq bam files first,
# as snakemake won't not know beforehand which ChIP-seq bam files belong to the TF x cell type combo.
rule count_normalize_chip:
    input:
        sites = SITES_FILE,
        chrom_size = CHROM_SIZE_FILE,
        chip = ALL_CHIP_BAM_FILES,
        chip_idxstats = ALL_CHIP_IDXSTATS_FILES,
        metadata = METADATA_TABLE_FILE
    output:
        touch(CHIP_COUNTS_FILE)
    params:
        ref_size = CHIP_REF_SIZE
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/count_normalize_chip/{VER_GENOME}/{{cell_type}}/{{tf_name}}_{{pwm_id}}_{THRESH_PVALUE}.log'
    shell:
        '''
        module load R/3.6.0
        Rscript ../scripts/count_chipseq_coverage.R \
        --pwm_id {wildcards.pwm_id} \
        --cell_type {wildcards.cell_type} \
        --sites {input.sites} \
        --chip_dir {CHIP_BAM_DIR} \
        --chrom_size {input.chrom_size} \
        --metadata {input.metadata} \
        --ref_size {params.ref_size} \
        --outdir {CHIP_COUNTS_DIR} \
        --outname {wildcards.tf_name}_{wildcards.pwm_id}_{THRESH_PVALUE}.{wildcards.cell_type} \
        &> {log}
        '''

# Count ChIP-seq bigWig signals for candidate sites
rule count_chip_signals:
    input:
        sites = SITES_FILE,
        metadata = METADATA_TABLE_FILE
    output:
        touch(CHIP_SIGNALS_FILE)
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/count_chip_signals/{VER_GENOME}/{{cell_type}}/{{tf_name}}_{{pwm_id}}_{THRESH_PVALUE}.log'
    shell:
        '''
        module load R/3.6.0
        Rscript ../scripts/count_chipseq_bigwig_signals.R \
        --pwm_id {wildcards.pwm_id} \
        --cell_type {wildcards.cell_type} \
        --sites {input.sites} \
        --chip_dir {CHIP_BIGWIG_DIR} \
        --metadata {input.metadata} \
        --outdir {CHIP_SIGNALS_DIR} \
        --outname {wildcards.tf_name}_{wildcards.pwm_id}_{THRESH_PVALUE}.{wildcards.cell_type} \
        &> {log}
        '''

rule list_chip_samples:
    input: ALL_CHIP_BAM_FILES
    shell:
        ' echo "ChIP samples: {ALL_CHIP_SAMPLES}" '


rule all_chip_idxstats:
    input:
        ALL_CHIP_BAM_FILES,
        ALL_CHIP_BAI_FILES,
        ALL_CHIP_IDXSTATS_FILES


rule all_chip_counts:
    input:
        ALL_CHIP_BAM_FILES,
        ALL_CHIP_BAI_FILES,
        ALL_CHIP_IDXSTATS_FILES,
        ALL_CHIP_COUNTS_FILES
    shell:
        ' echo "ChIP samples: {ALL_CHIP_SAMPLES}" '

# rule all_chip_signals:
#     input:
#         ALL_CHIP_SIGNALS_FILES
#     shell:
#         ' echo "ChIP samples: {ALL_CHIP_SAMPLES}" '
