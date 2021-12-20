
# Count DNase cuts along the genome using bedtools and convert to bigWig (.bw) format
rule count_dnase_genome_cuts:
    input:
        bam = DNASE_BAM_FILE,
        chrom_size = CHROM_SIZE_FILE,
    output:
        fwd = DNASE_GENOMECOUNTS_FWD_FILE,
        rev = DNASE_GENOMECOUNTS_REV_FILE
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/count_dnase_genome_cuts/{VER_GENOME}/{{dnase_sample}}.log'
    shell:
        '''
        module load R/3.6.0
        Rscript ../scripts/count_dnase_genome_cuts.R \
        --bam {input.bam} \
        --chrom_size {input.chrom_size} \
        --out_dir {DNASE_GENOMECOUNTS_DIR} \
        --out_prefix {wildcards.dnase_sample} \
        &> {log}
        '''


# Get DNase count matrices for a TF motif's candidate sites
rule get_dnase_motif_counts:
    input:
        sites = SITES_FILE,
        fwd = DNASE_GENOMECOUNTS_FWD_FILE,
        rev = DNASE_GENOMECOUNTS_REV_FILE
    output:
        touch(DNASE_COUNTMATRIX_FILE)
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/get_dnase_motif_counts/{VER_GENOME}/{{dnase_sample}}/{{pwm_id}}_{THRESH_PVALUE}.log'
    shell:
        '''
        module load R/3.6.0
        Rscript ../scripts/get_dnase_motif_counts.R \
        --sites {input.sites} \
        --dnase_fwd {input.fwd} \
        --dnase_rev {input.rev} \
        --outdir {DNASE_COUNTMATRIX_DIR}/{wildcards.dnase_sample} \
        --outname {wildcards.pwm_id}_{THRESH_PVALUE} \
        &> {log}
        '''


# Normalize and bin the DNase counts
rule normalize_bin_dnase:
    input:
        dnase = DNASE_COUNTMATRIX_FILE,
        idxstats = DNASE_IDXSTATS_FILE
    output:
        touch(DNASE_NORMALIZED_COUNTMATRIX_FILE),
        touch(DNASE_BINS_FILE)
    params:
        transform = 'asinh',
        ref_size = DNASE_REF_SIZE,
        bin = BIN,
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/normalize_bin_dnase/{VER_GENOME}/{{dnase_sample}}/{{pwm_id}}_{THRESH_PVALUE}.log'
    shell:
        '''
        module load R/3.6.0
        Rscript ../scripts/normalize_bin_dnase.R \
        --dnase_counts {input.dnase} \
        --idxstats {input.idxstats} \
        --bin {params.bin} \
        --ref_size {params.ref_size} \
        --transform {params.transform} \
        --outdir {DNASE_COUNTMATRIX_DIR}/{wildcards.dnase_sample} \
        --outname {wildcards.pwm_id}_{THRESH_PVALUE} \
        &> {log}
        '''


# get normalized_countmatrix.rds filenames, which will be checked for the input of combine_dnase_chip_data rule
def get_dnase_files(wildcards):
    metadata_table = pd.read_csv(METADATA_TABLE_FILE, sep='\t')
    metadata_table = metadata_table[(metadata_table["tf_name"] == wildcards.tf_name) &
                                    (metadata_table["cell_type"] == wildcards.cell_type)]
    sample_names = ';'.join(metadata_table["dnase_file"].tolist()).split(';')
    dnase_files = expand(DNASE_NORMALIZED_COUNTMATRIX_FILE, dnase_sample = sample_names, allow_missing=True)
    return dnase_files

# combine DNase and ChIP training data for a TF in a cell type
rule combine_dnase_chip_data:
    input:
        metadata = METADATA_TABLE_FILE,
        dnase = get_dnase_files,
        chip = CHIP_COUNTS_FILE
    output:
        touch(COMBINED_DNASE_CHIP_DATA_FILE)
    params:
        bin = BIN,
        transform = 'asinh',
        thresh_pValue = THRESH_PVALUE,
        chip_dir = CHIP_COUNTS_DIR,
        dnase_dir = DNASE_COUNTMATRIX_DIR,
        dnase_idxstats_dir = DNASE_BAM_DIR,
        dnase_ref_size = DNASE_REF_SIZE,
        outdir = COMBINED_DNASE_CHIP_DATA_DIR
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/combine_dnase_chip_data/{VER_GENOME}/{DATA_TYPE}/{MOTIF_SET}/{{tf_name}}_{{pwm_id}}_{THRESH_PVALUE}.{{cell_type}}.log'
    shell:
        '''
        module load R/3.6.0
        Rscript ../scripts/combine_dnase_chip_data.R \
        --tf_name {wildcards.tf_name} \
        --cell_type {wildcards.cell_type} \
        --metadata {input.metadata} \
        --thresh_pValue {params.thresh_pValue} \
        --chip_counts {input.chip} \
        --dnase_dir {params.dnase_dir} \
        --dnase_idxstats_dir {params.dnase_idxstats_dir} \
        --dnase_ref_size {params.dnase_ref_size} \
        --bin {params.bin} \
        --transform {params.transform} \
        --outdir {params.outdir} \
        --outname {wildcards.tf_name}_{wildcards.pwm_id}_{params.thresh_pValue}.{wildcards.cell_type} \
        &> {log}
        '''


# add ChIP peaks and signals to combined data for a TF in a cell type
rule add_chip_peaks_signals_to_combined_dnase_data:
    input:
        metadata = METADATA_TABLE_FILE,
        combined_data = COMBINED_DNASE_CHIP_DATA_FILE,
    output:
        touch(COMBINED_DNASE_CHIP_PEAKS_SIGNALS_DATA_FILE)
    params:
        bin = BIN,
        thresh_pValue = THRESH_PVALUE,
        combined_data_dir = COMBINED_DNASE_CHIP_DATA_DIR,
        chip_dir = f'{CHIP_DATA_DIR}/{VER_GENOME}'
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/add_chip_peaks_signals_to_combined_dnase_data/{VER_GENOME}/{DATA_TYPE}/{MOTIF_SET}/{{tf_name}}_{{pwm_id}}_{THRESH_PVALUE}.{{cell_type}}.log'
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


rule list_dnase_samples:
    input: ALL_DNASE_BAM_FILES
    shell:
        ' echo "DNase samples: {ALL_DNASE_SAMPLES}" '


rule all_dnase_idxstats:
    input:
        ALL_DNASE_BAM_FILES,
        ALL_DNASE_BAI_FILES,
        ALL_DNASE_IDXSTATS_FILES


rule all_dnase_genomecounts:
    input:
        ALL_DNASE_BAM_FILES,
        ALL_DNASE_IDXSTATS_FILES,
        ALL_DNASE_GENOMECOUNTS_FWD_FILES,
        ALL_DNASE_GENOMECOUNTS_REV_FILES


rule all_dnase_countmatrices:
    input:
        ALL_DNASE_COUNTMATRIX_FILES


rule all_dnase_bins:
    input:
        ALL_DNASE_IDXSTATS_FILES,
        ALL_DNASE_BINS_FILES,
        ALL_DNASE_NORMALIZED_COUNTMATRIX_FILES


rule all_combined_dnase_chip_data:
    input:
        ALL_COMBINED_DNASE_CHIP_DATA_FILES


rule all_combined_dnase_chip_peaks_signals_data:
    input:
        ALL_COMBINED_DNASE_CHIP_PEAKS_SIGNALS_DATA_FILES


rule clean_dnase:
    shell:
        '''
        rm -rf {OUT_DIR}/dnase_genome_counts/{VER_GENOME}
        rm -rf {OUT_DIR}/dnase_motif_counts/{VER_GENOME}/{DATA_TYPE}/{MOTIF_SET}
        rm -rf {OUT_DIR}/combine_dnase_chip_data/{VER_GENOME}/{DATA_TYPE}/{MOTIF_SET}
        '''
