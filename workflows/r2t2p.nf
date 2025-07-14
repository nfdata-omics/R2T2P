/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                 } from '../modules/nf-core/multiqc/main'
include { STAR_ALIGN              } from '../modules/nf-core/star/align/main'

include { PREPARE_REF             } from '../subworkflows/local/prepare_ref'
include { PREPARE_FASTQ           } from '../subworkflows/local/prepare_fastq'
include { BAM_SORT_STATS_SAMTOOLS } from '../subworkflows/nf-core/bam_sort_stats_samtools'

include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_r2t2p_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow R2T2P {

    take:

    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta       // value channel: path(fasta)
    ch_gtf         // value channel: path(gtf)
    ch_gff         // value channel: path(gff)
    ch_star_index  // value channel: path(star_index)

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // 
    // SUBWORKFLOW: Prepare reference genome
    //
    PREPARE_REF (
        ch_fasta,
        ch_gtf,
        ch_gff,
        ch_star_index
    )
    ch_versions = ch_versions.mix(PREPARE_REF.out.versions)

    //
    // SUBWORKFLOW: Prepare FastQ files
    //
    PREPARE_FASTQ (
        ch_samplesheet
    )
    ch_reads = PREPARE_FASTQ.out.reads
    ch_versions = ch_versions.mix(PREPARE_REF.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_FASTQ.out.multiqc_files)

    //
    // Map reads with STAR
    //
    STAR_ALIGN (
        ch_reads,
        PREPARE_REF.out.star_index.map { [ [:], it ] },
        PREPARE_REF.out.gtf.map { [ [:], it ] },
        false, // star_ignore_sjdbgtf
        "", // seq_platform
        "" // seq_center
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( STAR_ALIGN.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)
    ch_multiqc_files  = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.stats)
        .mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat)
        .mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats)



    orig_bam            = STAR_ALIGN.out.bam                                 // channel: [ val(meta), path(bam)            ]
    log_final           = STAR_ALIGN.out.log_final                           // channel: [ val(meta), path(log_final)      ]
    log_out             = STAR_ALIGN.out.log_out                             // channel: [ val(meta), path(log_out)        ]
    log_progress        = STAR_ALIGN.out.log_progress                        // channel: [ val(meta), path(log_progress)   ]
    bam_sorted          = STAR_ALIGN.out.bam_sorted                          // channel: [ val(meta), path(bam)            ]
    fastq               = STAR_ALIGN.out.fastq                               // channel: [ val(meta), path(fastq)          ]
    tab                 = STAR_ALIGN.out.tab                                 // channel: [ val(meta), path(tab)            ]
    orig_bam_transcript = STAR_ALIGN.out.bam_transcript                      // channel: [ val(meta), path(bam)            ]

    bam                 = BAM_SORT_STATS_SAMTOOLS.out.bam             // channel: [ val(meta), path(bam) ]
    bai                 = BAM_SORT_STATS_SAMTOOLS.out.bai             // channel: [ val(meta), path(bai) ]
    stats               = BAM_SORT_STATS_SAMTOOLS.out.stats           // channel: [ val(meta), path(stats) ]
    flagstat            = BAM_SORT_STATS_SAMTOOLS.out.flagstat        // channel: [ val(meta), path(flagstat) ]
    idxstats            = BAM_SORT_STATS_SAMTOOLS.out.idxstats        // channel: [ val(meta), path(idxstats) ]


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'r2t2p_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
