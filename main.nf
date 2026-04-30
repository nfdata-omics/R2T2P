#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nfdata-omics/r2t2p
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nfdata-omics/r2t2p
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { R2T2P  } from './workflows/r2t2p'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_r2t2p_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_r2t2p_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_r2t2p_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta            = getGenomeAttribute('fasta')
params.gff              = getGenomeAttribute('gff')
params.gtf              = getGenomeAttribute('gtf')
params.star_index       = getGenomeAttribute('star')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFDATAOMICS_R2T2P {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    // Define channels for reference files and other input files
    ch_fasta       = params.fasta      ? channel.value(file(params.fasta, checkIfExists: true))      : channel.empty()
    ch_gtf         = params.gtf        ? channel.value(file(params.gtf, checkIfExists: true))        : channel.empty()
    ch_gff         = params.gff        ? channel.value(file(params.gff, checkIfExists: true))        : channel.empty()
    ch_star_index  = params.star_index ? channel.value(file(params.star_index, checkIfExists: true)) : channel.empty()
    ch_user_gtf    = params.user_provided_annotation ? channel.value(file(params.user_provided_annotation, checkIfExists: true)) : channel.empty()
    ch_fragpipe_workflow = params.fragpipe_workflow ? channel.value(file(params.fragpipe_workflow, checkIfExists: true)) : channel.empty()
    ch_tools_folder = params.fragpipe_tools_folder ? channel.value(file(params.fragpipe_tools_folder, checkIfExists: true)) : channel.empty()
    ch_diann_folder = params.fragpipe_diann_folder ? channel.value(file(params.fragpipe_diann_folder, checkIfExists: true)) : channel.empty()
    ch_fragpipe_manifest = params.fragpipe_manifest ? channel.value(file(params.fragpipe_manifest, checkIfExists: true)) : channel.empty()
    ch_fragpipe_annotation = params.fragpipe_annotation ? channel.value(file(params.fragpipe_annotation, checkIfExists: true)) : channel.empty()

    //
    // WORKFLOW: Run pipeline
    //
    R2T2P (
        samplesheet,
        ch_fasta,
        ch_gtf,
        ch_gff,
        ch_star_index,
        ch_user_gtf,
        ch_fragpipe_workflow,
        ch_tools_folder,
        ch_diann_folder,
        ch_fragpipe_manifest,
        ch_fragpipe_annotation,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
    multiqc_report = R2T2P.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFDATAOMICS_R2T2P (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFDATAOMICS_R2T2P.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
