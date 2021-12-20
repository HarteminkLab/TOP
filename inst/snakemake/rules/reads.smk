
# Generate index and idxstats for sorted reads
rule samtools_index_idxstats:
    input:
        '{dir_bam}/{sample}.bam',
    output:
        bai = '{dir_bam}/{sample}.bam.bai',
        idxstats = '{dir_bam}/{sample}.bam.idxstats.txt'
    conda:
        '../envs/top.yaml'
    shell:
        '''
        samtools index {input}
        samtools idxstats {input} > {output.idxstats}
        '''

# Download ENCODE bam file
rule download_encode_bam_file:
    output:
        '{dir_bam}/{sample}.bam'
    shell:
        '''
        cd {wildcards.dir_bam}
        curl -O -L https://www.encodeproject.org/files/{wildcards.sample}/@@download/{wildcards.sample}.bam
        '''

# Download ENCODE bed file
rule download_encode_bed_file:
    output:
        '{dir_bed}/{sample}.bed.gz'
    shell:
        '''
        cd {wildcards.dir_bed}
        curl -O -L https://www.encodeproject.org/files/{wildcards.sample}/@@download/{wildcards.sample}.bed.gz
        '''

# Download ENCODE bigwig file
rule download_encode_bigwig_file:
    output:
        '{dir_bigwig}/{sample}.bigwig.gz'
    shell:
        '''
        cd {wildcards.dir_bigwig}
        curl -O -L https://www.encodeproject.org/files/{wildcards.sample}/@@download/{wildcards.sample}.bigwig.gz
        '''
