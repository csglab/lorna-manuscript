#!/bin/bash
set -eo pipefail

eval "$(conda shell.bash hook)"
conda activate /scratch/asabe/envs/PSIsigma

CWD='/scratch/asabe/projects/cassette-exons/psi-sigma'
cd ${CWD}/data

# for f in m84127_240426_234003_s3.*; do a=${f/m84127_240426_234003_s3./}; b=${a/.flnc.sorted/.Aligned.sortedByCoord.out}; c=${b/_LM2_/_LM2_LR_rep}; d=${c/_Par_/_Par_LR_rep}; mv $f $d; done
# for f in *_rep1.Aligned.sortedByCoord.out.bam; do prefix=${f/_rep1.Aligned.sortedByCoord.out.bam/}; ls ${prefix}_rep*.Aligned.sortedByCoord.out.bam > ${prefix}.group.txt; done

# ref_annot="gencode.v43.annotation.gtf"
# (grep "^#" ${ref_annot}; grep -v "^#" ${ref_annot} | sort -k1,1 -k4,4n) > ${ref_annot/.gtf/.sorted.gtf}

num_threads=32

# for group in *_Par*.group.txt
for group in *_Par_LR.group.txt

do
    prefix_Par=$(basename ${group} .group.txt)
    prefix_LM2=${prefix_Par/_Par/_LM2}
    prefix_output=${prefix_Par/_Par/}

    type=1
    nread=16
    [[ $prefix_output == *"_LR"* ]] && type=2 && nread=4

    echo "${prefix_output}: Type ${type}: Nread: ${nread}: ${prefix_Par} vs. ${prefix_LM2}"

    perl /scratch/asabe/projects/cassette-exons/psi-sigma/PSI-Sigma-2.3/dummyai.pl \
        --gtf gencode.v43.annotation.sorted.gtf \
        --name ${prefix_output} \
        --type ${type} \
        --nread ${nread} \
        --output results/${prefix_output} \
        --groupa ${prefix_Par}.group.txt \
        --groupb ${prefix_LM2}.group.txt \
        --threads ${num_threads}
done
