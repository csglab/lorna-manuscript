process BamToFastq {
    label 'basic'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/flnc", mode: 'copy'

    input:
    path flnc_bam

    output:
    path flnc_fastq

    script:
    flnc_fastq = flnc_bam.getBaseName() + '.fastq'
    """
    samtools fastq ${flnc_bam} > ${flnc_fastq}
    """
}

process AlignFLNCReads {
    label 'multi_parallel'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/aligned", mode: 'copy'

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
        -uf \
        --MD \
        ${ref_genome_fasta} \
        ${flnc_fastq} \
        > ${flnc_aligned_sam}
    """
}

process GetSpliceJunctionsFromReferenceAnnotation {
    label 'basic'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/splice-junctions", mode: 'copy'

    input:
    path ref_annotation_gtf
    path ref_genome_fasta

    output:
    path ref_splice_junctions_tsv, emit: ref_splice_junctions_tsv

    script:
    ref_splice_junctions_tsv = ref_annotation_gtf.getBaseName() + '.splice-junctions.tsv'

    """
    mkdir -p $workDir/tmp
    export TMPDIR=$workDir/tmp
    python3.7 \$(which get_SJs_from_gtf.py) \
        --f ${ref_annotation_gtf} \
        --g ${ref_genome_fasta} \
        --o ${ref_splice_junctions_tsv}
    """
}

process CleanTranscripts {
    label 'multi_parallel'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/cleaned", mode: 'copy'

    input:
    path flnc_aligned_sam
    path ref_genome_fasta
    path ref_splice_junctions_tsv
    path ref_variation_vcf


    output:
    path flnc_cleaned_sam, emit: flnc_cleaned_sam
    path "${prefix}_clean.{fa,log,TE.log}", emit: other

    script:
    num_threads = task.cpus
    prefix = flnc_aligned_sam.getBaseName()
    flnc_cleaned_sam = prefix + '_clean.sam'

    """
    mkdir -p $workDir/tmp
    export TMPDIR=$workDir/tmp ## setting TMPDIR to for pybedtools temp files

    python $launchDir/scripts/TranscriptClean/TranscriptClean.py \
        --sam ${flnc_aligned_sam} \
        --genome ${ref_genome_fasta} \
        --threads ${num_threads} \
        --spliceJns ${ref_splice_junctions_tsv} \
        --correctMismatches true \
        --correctIndels true \
        --variants ${ref_variation_vcf} \
        --maxLenIndel 5 \
        --maxSJOffset 5 \
        --canonOnly \
        --tmpDir $workDir/tmp \
        --deleteTmp \
        --outprefix ${prefix}
    """
}

process LabelReadsForAsFraction {
    label 'multi_parallel'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/labelled", mode: 'copy'

    input:
    path flnc_cleaned_sam
    path ref_genome_fasta
    path ref_genome_fasta_index

    output:
    path flnc_labelled_sam, emit: flnc_labelled_sam
    path "*_read_labels.tsv", emit: read_labels_tsv

    script:
    num_threads = task.cpus
    prefix = flnc_cleaned_sam.getBaseName()
    flnc_labelled_sam = prefix + '_labeled.sam'

    """
    mkdir -p $workDir/tmp
    export TMPDIR=$workDir/tmp

    talon_label_reads \
        --f ${flnc_cleaned_sam} \
        --g ${ref_genome_fasta} \
        --t ${num_threads} \
        --ar 20 \
        --deleteTmp \
        --o ${prefix}
    """
}

process SetupTALONDatabase {
    label 'parallel'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/bulk_sample/talon-db", mode: 'copy'

    input:
    path ref_annotation_gtf
    path flnc_labelled_sams
    path bulk_sample_talon_config_csv

    output:
    path "bulk_sample.talon.db", emit: talon_db
    path "bulk_sample_QC.log", emit: qc_log
    path bulk_sample_talon_config_csv, emit: config_csv

    script:
    // 'gencode.v43.annotation.gtf'
    gencode_version = ref_annotation_gtf.getBaseName().split('[.]')[1]
    num_threads = task.cpus

    """
    mkdir -p $workDir/tmp
    export TMPDIR=$workDir/tmp

    talon_initialize_database \
        --f ${ref_annotation_gtf} \
        --a gencode_${gencode_version} \
        --g hg38 \
        --l 0 \
        --idprefix NOVEL \
        --5p 500 \
        --3p 300 \
        --o bulk_sample.talon

    talon \
        --f ${bulk_sample_talon_config_csv} \
        --db bulk_sample.talon.db \
        --build hg38 \
        --threads ${num_threads} \
        --cov 0.9 \
        --identity 0.8 \
        --o bulk_sample
    """
}

process QuantifyIsoformAbundance {
    label 'semi_basic'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/bulk_sample/isoform-abundance", mode: 'copy'

    input:
    path talon_db
    path ref_annotation_gtf

    output:
    path "*.tsv", emit: abundance_tsv

    script:
    gencode_version = ref_annotation_gtf.getBaseName().split('[.]')[1]

    """
    talon_abundance \
        --db ${talon_db} \
        -a gencode_${gencode_version} \
        --build hg38 \
        --o bulk_sample
    """

}

process QuantifyIsoformAbundanceWithFiltering {
    label 'semi_basic'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/bulk_sample/isoform-abundance-filtered", mode: 'copy'

    input:
    path talon_db
    path ref_annotation_gtf

    output:
    path "bulk_sample.whitelist.csv", emit: whitelist_csv
    path "*.tsv", emit: filtered_abundance_tsv

    script:
    gencode_version = ref_annotation_gtf.getBaseName().split('[.]')[1]

    """
    talon_filter_transcripts \
        --db ${talon_db} \
        -a gencode_${gencode_version} \
        --maxFracA 0.5 \
        --minCount 5 \
        --minDatasets 2 \
        --o bulk_sample.whitelist.csv

    talon_abundance \
        --db ${talon_db} \
        -a gencode_${gencode_version} \
        --build hg38 \
        --whitelist bulk_sample.whitelist.csv \
        --o bulk_sample
    """

}

process GenerateIsoformsWhiteList {
    label 'basic'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/bulk_sample/isoform-whitelist", mode: 'copy'

    input:
    val cell_line
    path talon_db
    path ref_annotation_gtf

    output:
    path "${cell_line}.whitelist.csv", emit: whitelist_csv

    script:
    gencode_version = ref_annotation_gtf.getBaseName().split('[.]')[1]

    """
    talon_filter_transcripts \
        --db ${talon_db} \
        -a gencode_${gencode_version} \
        --maxFracA 0.5 \
        --minCount 5 \
        --minDatasets 2 \
        --datasets ${cell_line}_1,${cell_line}_2 \
        --o ${cell_line}.whitelist.csv
    """

}

process QuantifyIsoformAbundanceWithFilteringPerCellLine {
    label 'basic'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/bulk_sample/isoform-abundance-filtered-per-cell-line", mode: 'copy'

    input:
    path whitelist_csvs
    path talon_db
    path ref_annotation_gtf

    output:
    path "*.tsv", emit: filtered_abundance_tsv
    path "all_cell_lines.whitelist.csv", emit: all_cell_lines_whitelist_csv

    script:
    gencode_version = ref_annotation_gtf.getBaseName().split('[.]')[1]

    """
    cat ${whitelist_csvs} | uniq > all_cell_lines.whitelist.csv

    talon_abundance \
        --db ${talon_db} \
        -a gencode_${gencode_version} \
        --build hg38 \
        --whitelist all_cell_lines.whitelist.csv \
        --o all_cell_lines
    """

}


process GenerateTranscriptomeAnnotation {
    label 'basic'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/bulk_sample/transcriptome-gtf", mode: 'copy'

    input:
    path talon_db
    path whitelist_csv
    path ref_annotation_gtf

    output:
    path "bulk_sample_talon.gtf", emit: transcriptome_gtf

    script:
    gencode_version = ref_annotation_gtf.getBaseName().split('[.]')[1]

    """
    talon_create_GTF \
        --db ${talon_db} \
        --annot gencode_${gencode_version} \
        --build hg38 \
        --whitelist ${whitelist_csv} \
        --o bulk_sample
    """

}

process ExtractTranscriptomeSequences {
    label 'basic'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/bulk_sample/transcriptome-fasta", mode: 'copy'

    input:
    path transcriptome_gtf
    path ref_genome_fasta
    path ref_genome_fasta_index

    output:
    path "bulk_sample.fasta", emit: transcriptome_fasta

    script:
    """
    gffread --gtf ${transcriptome_gtf} -g ${ref_genome_fasta} -w bulk_sample.fasta
    """
}


workflow {

    params.data_dir = params.project_dir + '/data/long-reads'
    params.ref_genome_fasta = params.project_dir + '/data/references/genome/GRCh38_no_alt_analysis_set_GCA_000001405.15.id_trimmed.fasta'
    params.ref_genome_fasta_index = params.ref_genome_fasta + '.fai'
    params.ref_annotation_gtf = params.project_dir + '/data/references/annotation/gencode.v43.annotation.gtf'

    // flnc_fastq_pattern = params.data_dir + '/raw/flnc/' + 'ENC*.flnc.sample.fastq'
    // flnc_fastqs = Channel
    //                 .fromPath(flnc_fastq_pattern)
    //                 .toSortedList()
    //                 .flatMap()

    // AlignFLNCReads(
    //     flnc_fastqs,
    //     params.ref_genome_fasta
    // )

    // ref_variation_vcf = params.project_dir + '/data/references/variation/00-common_all.vcf.gz'

    // GetSpliceJunctionsFromReferenceAnnotation(
    //     params.ref_annotation_gtf,
    //     params.ref_genome_fasta)

    // ref_splice_junctions_tsv = GetSpliceJunctionsFromReferenceAnnotation.out.ref_splice_junctions_tsv

    // CleanTranscripts(
    //     AlignFLNCReads.out.flnc_aligned_sam,
    //     params.ref_genome_fasta,
    //     ref_splice_junctions_tsv,
    //     ref_variation_vcf
    // )

    flnc_cleaned_sam_pattern = params.data_dir + '/batches-all/encode-pipeline/cleaned/' + '*.sample_clean.sam'
    flnc_cleaned_sams = Channel
                    .fromPath(flnc_cleaned_sam_pattern)
                    .toSortedList()
                    .flatten()
    flnc_cleaned_sams.view()

    LabelReadsForAsFraction(
    //     CleanTranscripts.out.flnc_cleaned_sam,
        flnc_cleaned_sams,
        params.ref_genome_fasta,
        params.ref_genome_fasta_index
    )

    // // flnc_labelled_sams
    LabelReadsForAsFraction.out.flnc_labelled_sam
        .map { sam_file ->
            cell_line = sam_file.getBaseName().split('[.]')[0]
            sam_filename = sam_file.getName()
            config_line = "${cell_line},${cell_line},PacBio-Sequel2,${sam_filename}"
            return config_line
        }
        .collectFile(name: 'talon_config.csv', newLine: true, sort: true)
        .set { bulk_sample_talon_config_csv }

    SetupTALONDatabase(
        params.ref_annotation_gtf,
        LabelReadsForAsFraction.out.flnc_labelled_sam.collect(),
        // flnc_labelled_sams.collect(),
        bulk_sample_talon_config_csv
    )

    QuantifyIsoformAbundance(
        SetupTALONDatabase.out.talon_db,
        params.ref_annotation_gtf
    )

    QuantifyIsoformAbundanceWithFiltering(
        SetupTALONDatabase.out.talon_db,
        params.ref_annotation_gtf
    )

    // bulk_sample_talon_db = params.data_dir + '/batches-all/encode-pipeline/bulk_sample/talon-db/bulk_sample.talon.db'
    // bulk_sample_filtered_abundance_tsv = params.data_dir + '/batches-all/encode-pipeline/bulk_sample/isoform-abundance-filtered/bulk_sample_talon_abundance_filtered.tsv'
    // bulk_sample_whitelist_csv = params.data_dir + '/batches-all/encode-pipeline/bulk_sample/isoform-abundance-filtered/bulk_sample.whitelist.csv'

    GenerateTranscriptomeAnnotation(
        SetupTALONDatabase.out.talon_db,
        QuantifyIsoformAbundanceWithFiltering.out.whitelist_csv,
        params.ref_annotation_gtf
    )

    ExtractTranscriptomeSequences(
        GenerateTranscriptomeAnnotation.out.transcriptome_gtf,
        params.ref_genome_fasta,
        params.ref_genome_fasta_index
    )

}