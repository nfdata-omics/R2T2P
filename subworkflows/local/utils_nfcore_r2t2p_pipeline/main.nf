//
// Subworkflow with functionality specific to the nfdata-omics/r2t2p pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ [meta.id, meta.assay_type], meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ [meta.id, meta.assay_type], meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()

    // either gff or gtf should be provided, not both
    if (!params.gff && !params.gtf) {
        error("Please check input parameters -> Neither GFF nor GTF files were provided. Please provide one of them.")
    }
    if (params.gff && params.gtf) {
        error("Please check input parameters -> Both GFF and GTF files were provided. Please provide only one of them.")
    }

    // when the fragpipe manifest is provided, all the additional fragpipe files must also be provided
    if (params.fragpipe_manifest && (!params.fragpipe_workflow || !params.fragpipe_TMT_annotation)) {
        error("Please check input parameters -> When fragpipe_manifest is provided, fragpipe_workflow and fragpipe_TMT_annotation must also be provided.")
    }

    // check that the required jars and executables for fragpipe are present in the folders provided
    if (params.fragpipe_manifest) {
        def ionquant_jar = file("${params.fragpipe_tools_folder}/IonQuant-*/IonQuant-*.jar")
        if (ionquant_jar.isEmpty()) {
            error("Please check input parameters -> IonQuant JAR file not found in ${params.fragpipe_tools_folder}/IonQuant-*/IonQuant-*.jar")
        }
        def msfragger_jar = file("${params.fragpipe_tools_folder}/MSFragger-*/MSFragger-*.jar")
        if (msfragger_jar.isEmpty()) {
            error("Please check input parameters -> MSFragger JAR file not found in ${params.fragpipe_tools_folder}/MSFragger-*/MSFragger-*.jar")
        }
        def diatracer_jar = file("${params.fragpipe_tools_folder}/diaTracer-*/diaTracer-*.jar")
        if (diatracer_jar.isEmpty()) {
            error("Please check input parameters -> diaTracer JAR file not found in ${params.fragpipe_tools_folder}/diaTracer-*/diaTracer-*.jar")
        }
        def diann_exe = file("${params.fragpipe_diann_folder}/diann-linux")
        if (!diann_exe.exists()) {
            error("Please check input parameters -> DIANN executable not found in ${params.fragpipe_diann_folder}/diann-linux")
        }
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "fq,",
            "FastQC (Andrews 2010),",
            "STAR (Dobin et al. 2013),",
            "SAMtools (Li et al. 2009),",
            "StringTie (Pertea et al. 2015),",
            "GffRead and GffCompare (Pertea and Pertea 2020),",
            "UCSC Genome Browser utilities (Lee et al. 2022),",
            "Bioconductor genomic infrastructure (Lawrence et al. 2013),",
            "Ribo-seQC/RiboseQC (Calviello et al. 2019),",
            "DESeq2 (Love et al. 2014),",
            "DEXSeq (Anders et al. 2012),",
            "ORFquant (Calviello et al. 2020; Harnett et al. 2021),",
            "MSFragger (Kong et al. 2017),",
            "Philosopher (da Veiga Leprevost et al. 2020),",
            "IonQuant (Yu et al. 2021),",
            "Percolator (Käll et al. 2007),",
            "DIA-NN (Demichev et al. 2020),",
            "FragPipe DIA and TMT workflows (Yu et al. 2023; Chang et al. 2026),",
            "and MultiQC (Ewels et al. 2016)."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>St. Jude Rust Labs. fq: FASTQ parser and linter. URL: <a href=\"https://github.com/stjude-rust-labs/fq\">https://github.com/stjude-rust-labs/fq</a></li>",
            "<li>Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. URL: <a href=\"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/\">https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</a></li>",
            "<li>Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15-21. doi: <a href=\"https://doi.org/10.1093/bioinformatics/bts635\">10.1093/bioinformatics/bts635</a></li>",
            "<li>Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078-2079. doi: <a href=\"https://doi.org/10.1093/bioinformatics/btp352\">10.1093/bioinformatics/btp352</a></li>",
            "<li>Pertea, M., Pertea, G. M., Antonescu, C. M., Chang, T. C., Mendell, J. T., & Salzberg, S. L. (2015). StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nature Biotechnology, 33(3), 290-295. doi: <a href=\"https://doi.org/10.1038/nbt.3122\">10.1038/nbt.3122</a></li>",
            "<li>Pertea, G., & Pertea, M. (2020). GFF Utilities: GffRead and GffCompare. F1000Research, 9, ISCB Comm J-304. doi: <a href=\"https://doi.org/10.12688/f1000research.23297.2\">10.12688/f1000research.23297.2</a></li>",
            "<li>Lee, B. T., Barber, G. P., Benet-Pagès, A., Casper, J., Clawson, H., Diekhans, M., Fischer, C., Navarro Gonzalez, J., Hinrichs, A. S., Lee, C. M., Muthuraman, P., Nassar, L. R., Nguy, B., Pereira, T., Perez, G., Raney, B. J., Rosenbloom, K. R., Schmelter, D., Speir, M. L., Wick, B., Zweig, A. S., Haussler, D., Kuhn, R. M., Haeussler, M., & Kent, W. J. (2022). The UCSC Genome Browser database: 2022 update. Nucleic Acids Research, 50(D1), D1115-D1122. doi: <a href=\"https://doi.org/10.1093/nar/gkab959\">10.1093/nar/gkab959</a></li>",
            "<li>Lawrence, M., Huber, W., Pagès, H., Aboyoun, P., Carlson, M., Gentleman, R., Morgan, M. T., & Carey, V. J. (2013). Software for Computing and Annotating Genomic Ranges. PLoS Computational Biology, 9(8), e1003118. doi: <a href=\"https://doi.org/10.1371/journal.pcbi.1003118\">10.1371/journal.pcbi.1003118</a></li>",
            "<li>Calviello, L., Sydow, D., Harnett, D., & Ohler, U. (2019). Ribo-seQC: comprehensive analysis of cytoplasmic and organellar ribosome profiling data. bioRxiv. doi: <a href=\"https://doi.org/10.1101/601468\">10.1101/601468</a></li>",
            "<li>Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550. doi: <a href=\"https://doi.org/10.1186/s13059-014-0550-8\">10.1186/s13059-014-0550-8</a></li>",
            "<li>Anders, S., Reyes, A., & Huber, W. (2012). Detecting differential usage of exons from RNA-seq data. Genome Research, 22(10), 2008-2017. doi: <a href=\"https://doi.org/10.1101/gr.133744.111\">10.1101/gr.133744.111</a></li>",
            "<li>Calviello, L., Hirsekorn, A., & Ohler, U. (2020). Quantification of translation uncovers the functions of the alternative transcriptome. Nature Structural & Molecular Biology, 27(8), 717-725. doi: <a href=\"https://doi.org/10.1038/s41594-020-0450-4\">10.1038/s41594-020-0450-4</a></li>",
            "<li>Harnett, D., Meerdink, E., Calviello, L., Sydow, D., & Ohler, U. (2021). Genome-Wide Analysis of Actively Translated Open Reading Frames Using RiboTaper/ORFquant. Methods in Molecular Biology, 2252, 331-346. doi: <a href=\"https://doi.org/10.1007/978-1-0716-1150-0_16\">10.1007/978-1-0716-1150-0_16</a></li>",
            "<li>Kong, A. T., Leprevost, F. V., Avtonomov, D. M., Mellacheruvu, D., & Nesvizhskii, A. I. (2017). MSFragger: ultrafast and comprehensive peptide identification in mass spectrometry-based proteomics. Nature Methods, 14(5), 513-520. doi: <a href=\"https://doi.org/10.1038/nmeth.4256\">10.1038/nmeth.4256</a></li>",
            "<li>da Veiga Leprevost, F., Haynes, S. E., Avtonomov, D. M., Chang, H. Y., Shanmugam, A. K., Mellacheruvu, D., Kong, A. T., & Nesvizhskii, A. I. (2020). Philosopher: a versatile toolkit for shotgun proteomics data analysis. Nature Methods, 17(9), 869-870. doi: <a href=\"https://doi.org/10.1038/s41592-020-0912-y\">10.1038/s41592-020-0912-y</a></li>",
            "<li>Yu, F., Haynes, S. E., & Nesvizhskii, A. I. (2021). IonQuant Enables Accurate and Sensitive Label-Free Quantification With FDR-Controlled Match-Between-Runs. Molecular & Cellular Proteomics, 20, 100077. doi: <a href=\"https://doi.org/10.1016/j.mcpro.2021.100077\">10.1016/j.mcpro.2021.100077</a></li>",
            "<li>Käll, L., Canterbury, J. D., Weston, J., Noble, W. S., & MacCoss, M. J. (2007). Semi-supervised learning for peptide identification from shotgun proteomics datasets. Nature Methods, 4(11), 923-925. doi: <a href=\"https://doi.org/10.1038/nmeth1113\">10.1038/nmeth1113</a></li>",
            "<li>Demichev, V., Messner, C. B., Vernardis, S. I., Lilley, K. S., & Ralser, M. (2020). DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. Nature Methods, 17(1), 41-44. doi: <a href=\"https://doi.org/10.1038/s41592-019-0638-x\">10.1038/s41592-019-0638-x</a></li>",
            "<li>Yu, F., Teo, G. C., Kong, A. T., Fröhlich, K., Li, G. X., Demichev, V., & Nesvizhskii, A. I. (2023). Analysis of DIA proteomics data using MSFragger-DIA and FragPipe computational platform. Nature Communications, 14, 4154. doi: <a href=\"https://doi.org/10.1038/s41467-023-39869-5\">10.1038/s41467-023-39869-5</a></li>",
            "<li>Chang, H. Y., Deng, Y., Li, R., Avtonomov, D., Wen, B., Haynes, S. E., da Veiga Leprevost, F., Zhang, B., Yu, F., & Nesvizhskii, A. I. (2026). Analysis of isobaric quantitative proteomic data using TMT-Integrator and FragPipe computational platform. Nature Communications. doi: <a href=\"https://doi.org/10.1038/s41467-026-70118-7\">10.1038/s41467-026-70118-7</a></li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
