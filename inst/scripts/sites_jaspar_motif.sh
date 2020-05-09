#!/bin/bash

################################################################################
# The script process FIMO motif matches to get candidate binding sites
#
# requires process_pwm_sites.R to get candidate sites from motif matches
# requires bigWigAverageOverBed tool from UCSC to compute mapablity
# requires bedtools to filter out ENCODE blacklist regions
################################################################################

## parameters
tf_name=$1
pwm_id=$2
thresh_pValue=$3 # 1e-4
dir_blacklist=$4
dir_output=$5

###### example ######
# tf_name="CTCF"
# pwm_id="MA0139.1"
# thresh_pValue=1e-4
#####################

ver_genome="hg19"

flank=100
max_fimo_sites=1000000

## extract pwm meme file
echo "pwm_id: $pwm_id, tf_name: $tf_name"

pwm_name="${tf_name}_${pwm_id}_${thresh_pValue}"

dir_code=$HOME/projects/ATAC-seq/ATAC-seq_workflow/code_RCC

# output directory for FIMO motif matching results
mkdir -p ${dir_output}

## FIMO matching motifs
dir_motif_matches="${dir_output}/FIMO/${pwm_name}"
filename_motif_matches="${dir_motif_matches}/fimo.txt"

## get candidate sites after taking flanking windows around motif matches

# output directory for candidate sites
dir_sites="${dir_output}/candidate_sites/${thresh_pValue}"
mkdir -p ${dir_sites}

sites="${dir_sites}/${pwm_name}_flank${flank}_fimo_sites.bed"

## mapability downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/
echo "Flank motif sites and compute mapability..."
filename_mapability="${dir_blacklist}/wgEncodeDukeMapabilityUniqueness20bp.bigWig"
Rscript process_pwm_sites.R ${tf_name} ${pwm_id} ${filename_motif_matches} ${flank} ${thresh_pValue} ${dir_sites} ${filename_mapability}

## filter ENCODE blacklist regions
echo "Filter ENCODE blacklist regions..."
blacklist_filename="${dir_blacklist}/wgEncodeDacMapabilityConsensusExcludable.bed"
bedtools intersect -v -a ${sites} -b ${blacklist_filename} > ${sites}.filtered.tmp
head -1 ${sites} | cat - ${sites}.filtered.tmp > ${sites}.filtered
mv ${sites}.filtered ${sites}
rm ${sites}.filtered.tmp

echo "Complete candidate sites processing."
