
# assemble training data
rule assemble_dnase_training_data:
    input:
        metadata = METADATA_TABLE_FILE,
        combined_data_files = ALL_COMBINED_DNASE_CHIP_DATA_FILES
    output:
        ASSEMBLED_TRAINING_DATA_FILES
    params:
        bin = BIN,
        transform = 'asinh',
        thresh_pValue = THRESH_PVALUE,
        outdir = TRAINING_DATA_DIR,
        tf_list = ','.join(TRAINING_TF_SET),
        celltype_list = ','.join(TRAINING_CELLTYPE_SET),
        combined_data_dir = COMBINED_DATA_DIR
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/assembled_training_data/{VER_GENOME}/{DATA_TYPE}/{MOTIF_SET}/assemble_training_data.log'
    shell:
        '''
        module load R/3.6.0
        Rscript ../scripts/assemble_training_data.R \
        --tf_list {params.tf_list} \
        --celltype_list {params.celltype_list} \
        --metadata {input.metadata} \
        --combined_data_dir {params.combined_data_dir} \
        --thresh_pValue {params.thresh_pValue} \
        --bin {params.bin} \
        --transform {params.transform} \
        --outdir {params.outdir} \
        --outname {DATA_TYPE}_{VER_GENOME}_training_data \
        &> {log}
        '''


rule all_combined_data:
    input:
        ALL_COMBINED_DATA_FILES

rule assembled_data:
    input:
        ASSEMBLED_TRAINING_DATA_FILES
