
# Get motif matches using FIMO
rule fimo:
    input:
        pwm = PWM_FILE,
        ref_genome = REF_GENOME_FILE
    output:
        touch(FIMO_FILE)
    params:
        max_stored_scores = 1000000,
        thresh_pValue = THRESH_PVALUE,
        outdir = FIMO_DIR
    log:
        f'{LOG_DIR}/fimo/{VER_GENOME}/{{pwm_id}}_{THRESH_PVALUE}.log'
    shell:
        '''
        mkdir -p {params.outdir}
        fimo --text \
        --skip-matched-sequence --verbosity 2 \
        --oc {params.outdir} \
        --bgfile --uniform-- \
        --thresh {params.thresh_pValue} \
        --max-stored-scores {params.max_stored_scores} \
        {input.pwm} {input.ref_genome} \
        > {output}
        '''


# Get candidate sites using FIMO motif matches
rule candidate_sites:
    input:
        fimo = FIMO_FILE,
        blacklist = BLACKLIST_FILE
    output:
        touch(SITES_FILE)
    params:
        flank = 100,
        thresh_pValue = THRESH_PVALUE
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/candidate_sites/{VER_GENOME}/{{pwm_id}}_{THRESH_PVALUE}.log'
    shell:
        '''
        module load R/3.6.0
        Rscript ../scripts/get_candidate_sites.R \
        --fimo {input.fimo} --blacklist {input.blacklist} \
        --flank {params.flank} --thresh_pValue {params.thresh_pValue} \
        --out {output} \
        &> {log}
        '''


rule list_pwms:
    input: ALL_PWM_FILES
    shell:
        ' echo "PWMs: {ALL_PWM_IDS}" '

rule all_fimo:
    input: ALL_FIMO_FILES

rule all_sites:
    input: ALL_SITES_FILES
