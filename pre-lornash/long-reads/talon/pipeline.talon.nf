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
    label 'parallel'
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
    label 'multi_parallel'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/talon-db", mode: 'copy'

    input:
    tuple val(cell_line), path(sam_files)
    path ref_annotation_gtf

    output:
    tuple val(cell_line), path("${cell_line}.talon.db"), emit: talon_db
    path "${cell_line}.config.csv", emit: config_csv
    path "${cell_line}_QC.log", emit: qc_log

    script:
    sam_files = sam_files.sort { it.name }
    rep1_sam = sam_files[0]
    rep2_sam = sam_files[1]

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
        --o ${cell_line}.talon

    touch ${cell_line}.config.csv
    echo "Rep1,${cell_line},PacBio-Sequel2,${rep1_sam}" >> ${cell_line}.config.csv
    echo "Rep2,${cell_line},PacBio-Sequel2,${rep2_sam}" >> ${cell_line}.config.csv

    talon \
        --f ${cell_line}.config.csv \
        --db ${cell_line}.talon.db \
        --build hg38 \
        --threads ${num_threads} \
        --cov 0.9 \
        --identity 0.8 \
        --o ${cell_line}

    """
}

process QuantifyIsoformAbundance {
    label 'basic'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/isoform-abundance", mode: 'copy'

    input:
    tuple val(cell_line), path(talon_db)
    path ref_annotation_gtf

    output:
    tuple val(cell_line), path("${cell_line}.whitelist.csv"), emit: whitelist_csv
    path "${cell_line}_talon_abundance_filtered.tsv", emit: abundance_tsv

    script:
    gencode_version = ref_annotation_gtf.getBaseName().split('[.]')[1]

    """
    talon_filter_transcripts \
        --db ${talon_db} \
        -a gencode_${gencode_version} \
        --maxFracA 0.5 \
        --minCount 5 \
        --minDatasets 2 \
        --datasets Rep1,Rep2 \
        --o ${cell_line}.whitelist.csv

    talon_abundance \
        --db ${talon_db} \
        -a gencode_${gencode_version} \
        --build hg38 \
        --whitelist ${cell_line}.whitelist.csv \
        --o ${cell_line}

        # -d Rep1,Rep2
    """

}

process GenerateTranscriptomeAnnotation {
    label 'basic'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/transcriptome-gtf", mode: 'copy'

    input:
    tuple val(cell_line), path(talon_db), path(whitelist_csv)
    path ref_annotation_gtf

    output:
    tuple val(cell_line), path("${cell_line}_talon.gtf"), emit: transcriptome_gtf

    script:
    gencode_version = ref_annotation_gtf.getBaseName().split('[.]')[1]

    """
    talon_create_GTF \
        --db ${talon_db} \
        --annot gencode_${gencode_version} \
        --build hg38 \
        --whitelist ${whitelist_csv} \
        --o ${cell_line}
    """

}

process ExtractTranscriptomeSequences {
    label 'basic'
    publishDir "${params.data_dir}/batches-all/encode-pipeline/transcriptome-fasta", mode: 'copy'

    input:
    tuple val(cell_line), path(transcriptome_gtf)
    path ref_genome_fasta
    path ref_genome_fasta_index

    output:
    tuple val(cell_line), path("${cell_line}.fasta"), emit: transcriptome_fasta

    script:
    """
    gffread --gtf ${transcriptome_gtf} -g ${ref_genome_fasta} -w ${cell_line}.fasta
    """
}

process TestProcess {

    input:
    val salam

    script:
    work_dir = task.workDir
    """
    echo $task > task
    echo ${task.workDir} > work_dir
    echo $work_dir > work_dir2
    echo $workDir > workDir3
    echo $projectDir > projectDir
    echo $launchDir > launchDir
    echo $moduleDir > moduleDir
    """

}


workflow {

    params.data_dir = params.project_dir + '/data/long-reads'
    params.ref_genome_fasta = params.project_dir + '/data/references/genome/GRCh38_no_alt_analysis_set_GCA_000001405.15.id_trimmed.fasta'
    params.ref_genome_fasta_index = params.ref_genome_fasta + '.fai'
    params.ref_annotation_gtf = params.project_dir + '/data/references/annotation/gencode.v43.annotation.gtf'
    
    // flnc_bam_pattern = params.data_dir + '/batch{1,2,3,4,5,6,7}/encode-pipeline/flnc/' + '*.flnc.bam'
    // flnc_bams = Channel
    //                 .fromPath(flnc_bam_pattern)
    //                 .toSortedList()
    //                 .flatMap()

    // BamToFastq(flnc_bams)

    flnc_fastq_pattern = params.data_dir + '/raw/flnc/' + '*.flnc.fastq'
    flnc_fastqs = Channel
                    .fromPath(flnc_fastq_pattern)
                    .toSortedList()
                    .flatMap()

    AlignFLNCReads(
        flnc_fastqs,
        params.ref_genome_fasta
    )

    ref_variation_vcf = params.project_dir + '/data/references/variation/00-common_all.vcf.gz'

    GetSpliceJunctionsFromReferenceAnnotation(
        params.ref_annotation_gtf,
        params.ref_genome_fasta)

    ref_splice_junctions_tsv = GetSpliceJunctionsFromReferenceAnnotation.out.ref_splice_junctions_tsv

    CleanTranscripts(
        AlignFLNCReads.out.flnc_aligned_sam,
        params.ref_genome_fasta,
        ref_splice_junctions_tsv,
        ref_variation_vcf
    )

    LabelReadsForAsFraction(
        CleanTranscripts.out.flnc_cleaned_sam,
        params.ref_genome_fasta,
        params.ref_genome_fasta_index
    )

    LabelReadsForAsFraction.out.flnc_labelled_sam
        .flatten()
        .map { sam_file ->
            def cell_line = sam_file.getBaseName().split('[.]')[1].replaceFirst(/_1$|_2$/, '')
            return tuple(cell_line, sam_file)
        }
        .groupTuple()
        .set { cell_line_pair_sams }

    SetupTALONDatabase(
        cell_line_pair_sams,
        params.ref_annotation_gtf
    )

    QuantifyIsoformAbundance(
        SetupTALONDatabase.out.talon_db,
        params.ref_annotation_gtf
    )

    SetupTALONDatabase.out.talon_db
        .combine(
            QuantifyIsoformAbundance.out.whitelist_csv,
            by: 0)
        .set { cell_line_db_whitelist_pairs } 

    GenerateTranscriptomeAnnotation(
        cell_line_db_whitelist_pairs,
        params.ref_annotation_gtf
    )

    ExtractTranscriptomeSequences(
        GenerateTranscriptomeAnnotation.out.transcriptome_gtf,
        params.ref_genome_fasta,
        params.ref_genome_fasta_index
    )

}