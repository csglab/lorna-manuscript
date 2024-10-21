#!/bin/bash
set -eo pipefail

ref_annot_gtf='references/gencode.v46.annotation.gtf'
cell_line='MDA231'
b1="samples/${cell_line}_Par.short_read.bams.txt"
b2="samples/${cell_line}_LM2.short_read.bams.txt"
od='post'
tmp='prep'
read_type='single'
read_length=76
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
