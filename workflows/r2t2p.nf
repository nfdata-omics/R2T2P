/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PREPARE_REF            } from '../subworkflows/local/prepare_ref'
include { PREPARE_FASTQ          } from '../subworkflows/local/prepare_fastq'
include { TWO_PASS_ALIGNMENT     } from '../subworkflows/local/two_pass_alignment/main'
include { TRANSCRIPTOME_ASSEMBLY } from '../subworkflows/local/transcriptome_assembly/main'
include { FINAL_ALIGNMENT        } from '../subworkflows/local/final_alignment/main'
include { DIFFERENTIAL_ANALYSIS  } from '../subworkflows/local/differential_analysis/main'
include { ORF_ANALYSIS           } from '../subworkflows/local/orf_analysis/main'
include { PROTEOMICS             } from '../subworkflows/local/proteomics/main'

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
    ch_star_index  // value channel: path(star_index)
    ch_user_gtf    // value channel: path(user_gtf)
    ch_fragpipe_workflow // value channel: path(fragpipe_workflow)
    ch_tools_folder // value channel: path(tools_folder)
    ch_diann_folder // value channel: path(diann_folder)
    ch_fragpipe_manifest // value channel: path(fragpipe_manifest)
    ch_fragpipe_annotation // value channel: path(fragpipe_TMT_annotation)
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    //
    // SUBWORKFLOW: Prepare reference genome
    //
    PREPARE_REF (
        ch_fasta,
        ch_gtf,
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
    ch_versions = ch_versions.mix(PREPARE_FASTQ.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_FASTQ.out.multiqc_files)

    ch_reads
        .branch { meta, _fastq ->
            rna:  meta.assay_type == "RNA"
            ribo: meta.assay_type == "Ribo"
        }
    .set { ch_reads_by_type }

    //
    // SUBWORKFLOW: Two-pass alignment with STAR
    //
    TWO_PASS_ALIGNMENT (
        ch_reads_by_type.rna,
        PREPARE_REF.out.star_index,
        PREPARE_REF.out.gtf,
        PREPARE_REF.out.fasta,
        PREPARE_REF.out.fai,
        PREPARE_REF.out.chrom_sizes,
        ch_samplesheet,
        PREPARE_REF.out.bsgenome,
        PREPARE_REF.out.gtf_Rannot
    )

    TWO_PASS_ALIGNMENT.out.strandedness
        .map { meta, strandedness_file ->

            def inferred_value = strandedness_file.readLines()
                .find { l -> l.startsWith('inferred_strands') }
                .split('\t')[1].trim()

            def strandedness = (inferred_value == '+') ? 'forward' :
                (inferred_value == '-') ? 'reverse' : 'unstranded'

            [meta.id, meta + [strandedness: strandedness]]
        }
        .join(
            TWO_PASS_ALIGNMENT.out.bam.map { meta, bam -> [meta.id, bam] }
        )
        .map { _id, updated_meta, bam -> [updated_meta, bam] }
        .set { ch_bam_with_strandness }

    //
    // Transcriptome assembly from aligned reads
    //
    TRANSCRIPTOME_ASSEMBLY (
        ch_bam_with_strandness,
        PREPARE_REF.out.fasta,
        PREPARE_REF.out.fai,
        PREPARE_REF.out.gtf,
        PREPARE_REF.out.bsgenome,
        PREPARE_REF.out.gtf_Rannot,
        ch_user_gtf
    )
    ch_versions = ch_versions.mix(TRANSCRIPTOME_ASSEMBLY.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(TRANSCRIPTOME_ASSEMBLY.out.gff_stats)

    ch_final_gtf = TRANSCRIPTOME_ASSEMBLY.out.gtf
        .map { _meta, file -> file }
        .toList()
        .map { [files: it] }
        .combine(
            PREPARE_REF.out.gtf
                .toList()
                .map { [files: it] }
        )
        .map { assembled, reference ->
            assembled.files ?: reference.files
        }
        .collect()

    ch_final_gtf_Rannot = TRANSCRIPTOME_ASSEMBLY.out.gtf_Rannot
        .toList()
        .map { [files: it] }
        .combine(
            PREPARE_REF.out.gtf_Rannot
                .toList()
                .map { [files: it] }
        )
        .map { assembled, reference ->
            (assembled.files ?: reference.files)[0]
        }
        .collect()

    //
    // Final alignment
    //
    FINAL_ALIGNMENT (
        ch_reads,
        PREPARE_REF.out.star_index,                                     // genome_index
        ch_final_gtf,                                                   // new gtf after assembly
        PREPARE_REF.out.fasta,                                          // genome fasta
        PREPARE_REF.out.fai,                                            // genome fai
        PREPARE_REF.out.chrom_sizes,                                    // chrom sizes for bigWig conversion
        PREPARE_REF.out.bsgenome,                                       // bsgenome for Ribo-seQC
        ch_final_gtf_Rannot                                             // gtf R-object for Ribo-seQC
    )
    ch_versions = ch_versions.mix(FINAL_ALIGNMENT.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(FINAL_ALIGNMENT.out.multiqc_files)

    //
    // Differential expression analysis
    //
    DIFFERENTIAL_ANALYSIS (
        FINAL_ALIGNMENT.out.counts_regions,                             // Ribo-seQC-processed alignment files
        PREPARE_REF.out.bsgenome,                                       // bsgenome for Ribo-seQC
        ch_final_gtf_Rannot                                             // gtf R-object for Ribo-seQC
    )
    ch_versions = ch_versions.mix(DIFFERENTIAL_ANALYSIS.out.versions)

    //
    // ORF analysis
    //
    ORF_ANALYSIS (
        FINAL_ALIGNMENT.out.bam,                                        // bam files from alignment against final transcriptome
        FINAL_ALIGNMENT.out.bam_for_orfquant,                           // RData files for ORFquant for the RiboSeq samples
        PREPARE_REF.out.bsgenome,                                       // bsgenome for Ribo-seQC
        ch_final_gtf_Rannot                                             // gtf R-object for Ribo-seQC
    )
    ch_versions = ch_versions.mix(ORF_ANALYSIS.out.versions)

    //
    // Proteomics database creation and analysis
    //
    PROTEOMICS (
        PREPARE_REF.out.bsgenome,
        PREPARE_REF.out.gtf_Rannot,
        ch_final_gtf_Rannot,
        ORF_ANALYSIS.out.orfquant_fasta,
        ch_fragpipe_workflow,
        ch_tools_folder,
        ch_diann_folder,
        ch_fragpipe_manifest,
        ch_fragpipe_annotation
    )
    ch_versions = ch_versions.mix(PROTEOMICS.out.versions)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'r2t2p_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'r2t2p'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
