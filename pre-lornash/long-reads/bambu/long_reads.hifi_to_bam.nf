process DemultiplexBarcodes {
    label 'parallel'

    input:
    path hifi_reads_bam
    path hifi_reads_bam_index
    path barcodes_fasta

    output:
    path "*.demux.*.bam",      emit: demux_bam
    path "*.demux.*.bam.pbi",  emit: demux_bam_index
    path "*.{xml,clips,json,counts,guess,report,summary}", emit: other

    script:
    num_threads = task.cpus
    // prefix = hifi_reads_bam.getBaseName().replace('.hifi_reads.default', '')
    prefix = hifi_reads_bam.getBaseName().replace('.hifi_reads', '')
    demux_bam = prefix + '.demux.bam'

    """
    lima \
        ${hifi_reads_bam} \
        ${barcodes_fasta} \
        ${demux_bam} \
        --ccs --same \
        --num-threads ${num_threads} \
        --split-named \
        --log-level INFO \
        --log-file ${prefix}.lima.demux.log.txt        
    """
        
}

process RemovePrimers {
    label 'multi_parallel'

    input:
    path demux_bam
    path demux_bam_index
    path primers_fasta

    output:
    path "*.rm_primers.*.bam",      emit: rm_primers_bam
    path "*.rm_primers.*.bam.pbi",  emit: rm_primers_bam_index
    path "*.{xml,clips,json,counts,guess,report,summary}", emit: other

    script:
    num_threads = task.cpus
    prefix = demux_bam.getBaseName()
    rm_primers_bam = prefix + '.rm_primers.bam'

    """
    lima \
        ${demux_bam} \
        ${primers_fasta} \
        ${rm_primers_bam} \
        --isoseq \
        --peek-guess \
        --num-threads ${num_threads} \
        --log-level INFO \
        --log-file ${prefix}.lima.rm_primers.log.txt
    """
}


process RenameSamples {
    label 'basic'

    input:
    path rm_primers_bam
    path rm_primers_bam_index

    output:
    path "*.fl.bam",      emit: fl_bam
    path "*.fl.bam.pbi",  emit: fl_bam_index

    script:
    num_threads = task.cpus
    // 'm64182_230108_073343.demux.A375tr_2--A375tr_2.rm_primers.NEB_5p--NEB_Clontech_3p'
    (pool_id, _, sample_id_twice, _, _) = rm_primers_bam.getBaseName().split('[.]')
    sample_id = sample_id_twice.split('--')[0]

    fl_prefix = [pool_id, sample_id, 'fl'].join('.')
    fl_bam = fl_prefix + '.bam'
    fl_bam_index = fl_bam + '.pbi'

    """
    mv ${rm_primers_bam} ${fl_bam}
    mv ${rm_primers_bam_index} ${fl_bam_index}
    """

}

process RefineFullLengthReads {
    label 'multi_parallel'

    input:
    path fl_bam
    path fl_bam_index
    path primers_fasta

    output:
    path "*.flnc.bam",      emit: flnc_bam
    path "*.flnc.bam.pbi",  emit: flnc_bam_index
    path "*.{xml,json,csv}", emit: other

    script:
    num_threads = task.cpus
    prefix = fl_bam.getBaseName().replace('.fl', '')
    flnc_bam = prefix + '.flnc.bam'

    """
    isoseq3 refine \
        ${fl_bam} \
        ${primers_fasta} \
        ${flnc_bam} \
        --require-polya \
        --min-polya-length 20 \
        --num-threads ${num_threads} \
        --log-level INFO \
        --log-file ${prefix}.refine.log.txt
    """

}

process BamToFastq {
    label 'basic'

    input:
    path flnc_bam

    output:
    path flnc_fastq, emit: flnc_fastq

    script:
    flnc_fastq = flnc_bam.getBaseName() + '.fastq'
    """
    samtools fastq ${flnc_bam} > ${flnc_fastq}
    """
}

process AlignFLNCReads {
    label 'multi_parallel'

    input:
    path flnc_fastq
    path ref_genome_fasta

    output:
    path flnc_aligned_sam, emit: flnc_aligned_sam

    script:
    num_threads = task.cpus
    flnc_aligned_sam = flnc_fastq.getBaseName() + '.sam'

    """
    minimap2 \
        -t ${num_threads} \
        -a \
        -x splice:hq \
        ${ref_genome_fasta} \
        ${flnc_fastq} \
        > ${flnc_aligned_sam}
    """
}

process CompressSortIndexAlignedReads {
    label 'basic'

    input:
    path aligned_sam

    output:
    path "${prefix}.sorted.bam", emit: aligned_bam
    path "${prefix}.sorted.bam.bai", emit: aligned_bam_index

    script:
    prefix = aligned_sam.getBaseName()
    """
    samtools view -b ${aligned_sam} > ${prefix}.bam
    samtools sort ${prefix}.bam > ${prefix}.sorted.bam
    samtools index ${prefix}.sorted.bam
    """
}


workflow {

    // hifi_reads_bam = params.batch_dir + '/hifi_reads/' + params.pool_id + '.hifi_reads.default.bam'
    hifi_reads_bam = params.batch_dir + '/hifi_reads/' + params.pool_id + '.hifi_reads.bam'

    hifi_reads_bam_index = hifi_reads_bam + '.pbi'

    barcodes_fasta = params.batch_dir + '/barcodes/' + params.pool_id + '.barcodes.fasta'
    primers_fasta = params.batch_dir + '/barcodes/' + params.pool_id + '.primers.fasta'

    DemultiplexBarcodes(
        hifi_reads_bam,
        hifi_reads_bam_index,
        barcodes_fasta)

    RemovePrimers(
        DemultiplexBarcodes.out.demux_bam.flatten(),
        DemultiplexBarcodes.out.demux_bam_index.flatten(),
        primers_fasta)

    RenameSamples(
        RemovePrimers.out.rm_primers_bam,
        RemovePrimers.out.rm_primers_bam_index)

    RefineFullLengthReads(
        RenameSamples.out.fl_bam,
        RenameSamples.out.fl_bam_index,
        primers_fasta)

    BamToFastq(RefineFullLengthReads.out.flnc_bam)

    AlignFLNCReads(
        BamToFastq.out.flnc_fastq,
        params.ref_genome_fasta
    )

    CompressSortIndexAlignedReads(
        AlignFLNCReads.out.flnc_aligned_sam
    ) 

}
