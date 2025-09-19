/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                             } from '../modules/nf-core/multiqc/main'
include { STAR_ALIGN as STAR_FIRST_ALIGN      } from '../modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_WITH_NOVEL_JUNCT } from '../modules/nf-core/star/align/main'
include { CREATE_FIRSTPASS_JUNCTIONS          } from '../modules/local/create_firstpass_junctions/main'

include { PREPARE_REF                                      } from '../subworkflows/local/prepare_ref'
include { PREPARE_FASTQ                                    } from '../subworkflows/local/prepare_fastq'
include { BAM_SORT_STATS_SAMTOOLS as FIRST_BAM_SORT_STATS  } from '../subworkflows/nf-core/bam_sort_stats_samtools'
include { BAM_SORT_STATS_SAMTOOLS as SECOND_BAM_SORT_STATS } from '../subworkflows/nf-core/bam_sort_stats_samtools'

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
        ch_star_index,
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
    STAR_FIRST_ALIGN (
        ch_reads,
        PREPARE_REF.out.star_index.map { [ [:], it ] },
        PREPARE_REF.out.gtf.map { [ [:], it ] },
        "$projectDir/assets/NO_FILE", // empty arguments for additional_junctions
        false, // star_ignore_sjdbgtf
        "", // seq_platform
        "" // seq_center
    )
    ch_versions = ch_versions.mix(STAR_FIRST_ALIGN.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(STAR_FIRST_ALIGN.out.log_final.collect{it[1]})

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    FIRST_BAM_SORT_STATS ( STAR_FIRST_ALIGN.out.bam, PREPARE_REF.out.fasta.map { [ [:], it ] } )
    ch_versions = ch_versions.mix(FIRST_BAM_SORT_STATS.out.versions)
    ch_multiqc_files  = ch_multiqc_files.mix( FIRST_BAM_SORT_STATS.out.stats.collect{it[1]} )
        .mix( FIRST_BAM_SORT_STATS.out.flagstat.collect{it[1]} )
        .mix( FIRST_BAM_SORT_STATS.out.idxstats.collect{it[1]} )

    //
    // Get the table of novel junctions from the first STAR alignment
    //
    CREATE_FIRSTPASS_JUNCTIONS(
        STAR_FIRST_ALIGN.out.pass1_spl_juc_tab,
        PREPARE_REF.out.bsgenome,
        PREPARE_REF.out.gtf_Rannot
    )

    // Match the reads wit the corresponding junction table
    ch_reads_with_junctions = ch_reads.join(CREATE_FIRSTPASS_JUNCTIONS.out.pass1_junctions, by: 0)
    // Split into two channels:
    ch_reads_ordered = ch_reads_with_junctions.map { it[0..1] }
    ch_junctions_ordered = ch_reads_with_junctions.map { it[2] }

    //
    // Second STAR alignment using novel junctions
    //
    STAR_WITH_NOVEL_JUNCT (
        ch_reads_ordered,
        PREPARE_REF.out.star_index.map { [ [:], it ] },
        PREPARE_REF.out.gtf.map { [ [:], it ] },
        ch_junctions_ordered, // channel for additional junctions
        false, // star_ignore_sjdbgtf
        "", // seq_platform
        "" // seq_center
    )
    ch_versions = ch_versions.mix(STAR_WITH_NOVEL_JUNCT.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(STAR_WITH_NOVEL_JUNCT.out.log_final.collect{it[1]})

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    SECOND_BAM_SORT_STATS ( STAR_WITH_NOVEL_JUNCT.out.bam, PREPARE_REF.out.fasta.map { [ [:], it ] } )
    ch_versions = ch_versions.mix(SECOND_BAM_SORT_STATS.out.versions)
    ch_multiqc_files  = ch_multiqc_files.mix(SECOND_BAM_SORT_STATS.out.stats.collect{it[1]} )
        .mix(SECOND_BAM_SORT_STATS.out.flagstat.collect{it[1]} )
        .mix(SECOND_BAM_SORT_STATS.out.idxstats.collect{it[1]} )

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
