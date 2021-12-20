
# Count DNase cuts along the genome using bedtools and convert to bigWig (.bw) format
rule count_atac_genome_cuts:
    input:
        bam = ATAC_BAM_FILE,
        chrom_size = CHROM_SIZE_FILE,
    output:
        fwd = ATAC_GENOMECOUNTS_FWD_FILE,
        rev = ATAC_GENOMECOUNTS_REV_FILE
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/count_atac_genome_cuts/{VER_GENOME}/{{atac_sample}}.log'
    shell:
        '''
        module load R/3.6.0
        Rscript ../scripts/count_dnase_genome_cuts.R \
        --bam {input.bam} \
        --chrom_size {input.chrom_size} \
        --out_dir {ATAC_GENOMECOUNTS_DIR} \
        --out_prefix {wildcards.atac_sample} \
        &> {log}
        '''


# Get DNase count matrices for a TF motif's candidate sites
rule get_atac_motif_counts:
    input:
        sites = SITES_FILE,
        fwd = ATAC_GENOMECOUNTS_FWD_FILE,
        rev = ATAC_GENOMECOUNTS_REV_FILE
    output:
        touch(ATAC_COUNTMATRIX_FILE)
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/get_atac_motif_counts/{VER_GENOME}/{{atac_sample}}/{{pwm_id}}_{THRESH_PVALUE}.log'
    shell:
        '''
        module load R/3.6.0
        Rscript ../scripts/get_dnase_motif_counts.R \
        --sites {input.sites} \
        --dnase_fwd {input.fwd} \
        --dnase_rev {input.rev} \
        --outdir {ATAC_COUNTMATRIX_DIR}/{wildcards.atac_sample} \
        --outname {wildcards.pwm_id}_{THRESH_PVALUE} \
        &> {log}
        '''


# Normalize and bin the DNase counts
rule normalize_bin_atac:
    input:
        atac = ATAC_COUNTMATRIX_FILE,
        idxstats = ATAC_IDXSTATS_FILE
    output:
        touch(ATAC_NORMALIZED_COUNTMATRIX_FILE),
        touch(ATAC_BINS_FILE)
    params:
        transform = 'asinh',
        ref_size = ATAC_REF_SIZE,
        bin = BIN,
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/normalize_bin_atac/{VER_GENOME}/{{atac_sample}}/{{pwm_id}}_{THRESH_PVALUE}.log'
    shell:
        '''
        module load R/3.6.0
        Rscript ../scripts/normalize_bin_dnase.R \
        --dnase_counts {input.atac} \
        --idxstats {input.idxstats} \
        --bin {params.bin} \
        --ref_size {params.ref_size} \
        --transform {params.transform} \
        --outdir {ATAC_COUNTMATRIX_DIR}/{wildcards.atac_sample} \
        --outname {wildcards.pwm_id}_{THRESH_PVALUE} \
        &> {log}
        '''

# get normalized_countmatrix.rds filenames, which will be checked for the input of combine_atac_chip_data rule
def get_atac_files(wildcards):
    metadata_table = pd.read_csv(METADATA_TABLE_FILE, sep='\t')
    metadata_table = metadata_table[(metadata_table["tf_name"] == wildcards.tf_name) &
                                    (metadata_table["cell_type"] == wildcards.cell_type)]
    sample_names = ';'.join(metadata_table["atac_file"].tolist()).split(';')
    atac_files = expand(ATAC_NORMALIZED_COUNTMATRIX_FILE, atac_sample = sample_names, allow_missing=True)
    return atac_files

# combine ATAC and ChIP training data for a TF in a cell type
rule combine_atac_chip_data:
    input:
        metadata = METADATA_TABLE_FILE,
        atac = get_atac_files,
        chip = CHIP_COUNTS_FILE
    output:
        touch(COMBINED_ATAC_CHIP_DATA_FILE)
    params:
        bin = BIN,
        transform = 'asinh',
        thresh_pValue = THRESH_PVALUE,
        chip_dir = CHIP_COUNTS_DIR,
        atac_dir = ATAC_COUNTMATRIX_DIR,
        atac_idxstats_dir = ATAC_BAM_DIR,
        atac_ref_size = ATAC_REF_SIZE,
        outdir = COMBINED_ATAC_CHIP_DATA_DIR
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/combine_atac_chip_data/{VER_GENOME}/{DATA_TYPE}/{MOTIF_SET}/{{tf_name}}_{{pwm_id}}_{THRESH_PVALUE}.{{cell_type}}.log'
    shell:
        '''
        module load R/3.6.0
        Rscript ../scripts/combine_dnase_chip_data.R \
        --tf_name {wildcards.tf_name} \
        --cell_type {wildcards.cell_type} \
        --metadata {input.metadata} \
        --thresh_pValue {params.thresh_pValue} \
        --chip_counts {input.chip} \
        --dnase_dir {params.atac_dir} \
        --dnase_idxstats_dir {params.atac_idxstats_dir} \
        --dnase_ref_size {params.atac_ref_size} \
        --bin {params.bin} \
        --transform {params.transform} \
        --outdir {params.outdir} \
        --outname {wildcards.tf_name}_{wildcards.pwm_id}_{params.thresh_pValue}.{wildcards.cell_type} \
        &> {log}
        '''

# add ChIP peaks and signals to combined data for a TF in a cell type
rule add_chip_peaks_signals_to_combined_atac_data:
    input:
        metadata = METADATA_TABLE_FILE,
        combined_data = COMBINED_ATAC_CHIP_DATA_FILE,
    output:
        touch(COMBINED_ATAC_CHIP_PEAKS_SIGNALS_DATA_FILE)
    params:
        bin = BIN,
        thresh_pValue = THRESH_PVALUE,
        combined_data_dir = COMBINED_ATAC_CHIP_DATA_DIR,
        chip_dir = f'{CHIP_DATA_DIR}/{VER_GENOME}'
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/add_chip_peaks_signals_to_combined_atac_data/{VER_GENOME}/{DATA_TYPE}/{MOTIF_SET}/{{tf_name}}_{{pwm_id}}_{THRESH_PVALUE}.{{cell_type}}.log'
    shell:
        '''
        module load R/3.6.0
        Rscript ../scripts/add_chip_peaks_signals_to_combined_data.R \
        --tf_name {wildcards.tf_name} \
        --cell_type {wildcards.cell_type} \
        --metadata {input.metadata} \
        --combined_data_dir {params.combined_data_dir} \
        --chip_dir {params.chip_dir} \
        --thresh_pValue {params.thresh_pValue} \
        --bin {params.bin} \
        &> {log}
        '''

rule list_atac_samples:
    input: ALL_ATAC_BAM_FILES
    shell:
        ' echo "ATAC samples: {ALL_ATAC_SAMPLES}" '


rule all_atac_idxstats:
    input:
        ALL_ATAC_BAM_FILES,
        ALL_ATAC_BAI_FILES,
        ALL_ATAC_IDXSTATS_FILES


rule all_atac_genomecounts:
    input:
        ALL_ATAC_BAM_FILES,
        ALL_ATAC_IDXSTATS_FILES,
        ALL_ATAC_GENOMECOUNTS_FWD_FILES,
        ALL_ATAC_GENOMECOUNTS_REV_FILES


rule all_atac_countmatrices:
    input:
        ALL_ATAC_COUNTMATRIX_FILES


rule all_atac_bins:
    input:
        ALL_ATAC_IDXSTATS_FILES,
        ALL_ATAC_BINS_FILES,
        ALL_ATAC_NORMALIZED_COUNTMATRIX_FILES


rule all_combined_atac_chip_data:
    input:
        ALL_COMBINED_ATAC_CHIP_DATA_FILES


rule all_combined_atac_chip_peaks_signals_data:
    input:
        ALL_COMBINED_ATAC_CHIP_PEAKS_SIGNALS_DATA_FILES


rule clean_atac:
    shell:
        '''
        rm -rf {OUT_DIR}/atac_genome_counts/{VER_GENOME}
        rm -rf {OUT_DIR}/atac_motif_counts/{VER_GENOME}/{MOTIF_SET}
        rm -rf {OUT_DIR}/combine_atac_chip_data/{VER_GENOME}/{MOTIF_SET}
        '''
