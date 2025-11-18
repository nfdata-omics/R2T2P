process MULTI_DE {
    tag "${meta.id}"
    label 'process_low'

    container "docker.io/nfdata/riboseqc:v1.3.0-patched"

    input:
    tuple val(meta), path(counts_regions), path(table)
    path "R-user-lib/*"
    path gtf_Rannot
    val org_db_name
    path "R-user-lib/*"

    output:
    path "*.RData", emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def run_go = org_db_name ? 'T' : 'F'
    def go_package = org_db_name ? "gopckg = \"${org_db_name}\"" : ''
    def library_import = org_db_name ? "library(${org_db_name})" : ''
    """
    export R_LIBS_USER=\$PWD/R-user-lib

    Rscript --vanilla -e '
        library(RiboseQC)
        library(DESeq2)
        library(DEXSeq)
        library(topGO)
        library(randomForest)
        library(glmnet)
        library(ggrepel)
        library(GenomicFeatures)
        library(ORFik)
        ${library_import}

        RiboseQC:::prepare_multiDE(
            file_matrix = "${table}",
            annotation_file = "${gtf_Rannot}",
            dest_suffix = "${meta.id}",
            ${args}
        )

        RiboseQC:::run_multiDE(
            DEobj = "${meta.id}_multiDE_input.RData",
            annotation_file = "${gtf_Rannot}",
            cores = ${task.cpus},
            runGO= ${run_go},
            ${go_package}
            ${args2}
        )

        RiboseQC:::run_RF_regression(
            summary_multiDE = "${meta.id}_multiDE_results_summary.RData"
        )

        # Writing package versions to versions.yml
        x = sessionInfo()
        versions <- list(
            R = paste0(x\$R.version\$major, ".", x\$R.version\$minor)
        )
        for(i in seq_along(x\$otherPkgs)){
            pkg <- x\$otherPkgs[[i]]
            versions[[pkg\$Package]] <- pkg\$Version
        }
        versions <- list(MULTI_DE = versions)
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

    stub:
    """

    Rscript -e '
        library(RiboseQC)
        library(DESeq2)
        library(DEXSeq)

        # Writing package versions to versions.yml
        x = sessionInfo()
        versions <- list(
            R = paste0(x\$R.version\$major, ".", x\$R.version\$minor)
        )
        for(i in seq_along(x\$otherPkgs)){
            pkg <- x\$otherPkgs[[i]]
            versions[[pkg$Package]] <- pkg\$Version
        }
        versions <- list(MULTI_DE = versions)
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

}
