#!/bin/bash
set -eo pipefail

ref_annot_gtf='references/gencode.v46.annotation.gtf'
cell_line='HCC1806'
b1='samples/HCC1806_Par.short_read.bams.txt'
b2='samples/HCC1806_LM2.short_read.bams.txt'
od='post'
tmp='prep'
read_type='paired'
read_length=151
anchor_length=8
num_threads=128
task='both'

od="${od}/${cell_line}/short-read"
tmp="${tmp}/${cell_line}/short-read"
mkdir -p ${od} ${tmp}

rmats.py \
	--gtf ${ref_annot_gtf} \
	--b1 ${b1} \
	--b2 ${b2} \
	--od ${od} \
	--tmp ${tmp} \
	-t ${read_type} \
	--libType fr-unstranded \
	--readLength ${read_length} \
	--variable-read-length \
	--anchorLength ${anchor_length} \
	--nthread ${num_threads} \
	--task ${task} \
	--novelSS
