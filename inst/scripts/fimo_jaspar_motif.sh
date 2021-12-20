#!/bin/bash

##############################################################################
# Get motif matches using FIMO
#
# Requires FIMO from MEME suite
# FIMO instructions: http://meme-suite.org/doc/fimo.html?man_type=web
# MEME Suite installation: http://meme-suite.org/doc/install.html?man_type=web
##############################################################################

## parameters
tf_name=$1
pwm_id=$2
thresh_pValue=$3
ver_genome=$4
dir_pwm=$5
dir_ref_genome=$6
dir_output=$7

###### example ######
# tf_name="CTCF"
# pwm_id="MA0139.1"
# thresh_pValue="1e-4"
# ver_genome="hg19"
#####################

max_fimo_sites=1000000

## extract pwm meme file
echo "pwm_id: $pwm_id, tf_name: $tf_name"

pwm_name="${tf_name}_${pwm_id}_${thresh_pValue}"

# directory for Jaspar motifs
pwm_filename="${dir_pwm}/${pwm_id}.meme"

# reference genome fasta file
genome_filename="${dir_ref_genome}/${ver_genome}.fa"

## FIMO matching motifs
mkdir -p ${dir_output}
dir_motif_matches="${dir_output}/FIMO/${pwm_name}"
filename_motif_matches="${dir_motif_matches}/fimo.txt"

mkdir -p ${dir_motif_matches}
echo "Begin FIMO matching motif ..."

fimo --oc ${dir_motif_matches} --verbosity 2 --text \
    --bgfile --uniform-- --thresh ${thresh_pValue} \
    --max-stored-scores ${max_fimo_sites} \
    ${pwm_filename} ${genome_filename}  \
    > ${filename_motif_matches}

echo "Finish FIMO: ${dir_motif_matches}"
