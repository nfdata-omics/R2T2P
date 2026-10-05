process ORFDEX {
    tag "${meta.id}"
    label 'process_single', 'process_high_memory'

    container "docker.io/nfdata/riboseqc:v2.1.0-patched"

    input:
    tuple val(meta), path(bam_files), path(table)
    path orfquant_results

    output:
    path "*.RData", emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
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

        # RiboseQC 0.99.0 does not initialise the top-level ORFDEX result list
        # before assigning to ORFDEX[[i]]. Patch a local copy of the function
        # until the fix is available in the container image.
        run_ORFDEX <- RiboseQC:::run_ORFDEX
        run_ORFDEX_body <- as.list(body(run_ORFDEX))
        body(run_ORFDEX) <- as.call(c(
            run_ORFDEX_body[1],
            list(quote(ORFDEX <- list())),
            run_ORFDEX_body[-1]
        ))

        run_ORFDEX(
            ORFquant_res = "${orfquant_results}",
            file_matrix = "${table}",
            suffix = "${meta.id}",
            ${args}
            cores = ${task.cpus}
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
        versions <- list(ORFDEX = versions)
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
        library(topGO)
        library(randomForest)
        library(glmnet)
        library(ggrepel)
        library(GenomicFeatures)
        library(ORFik)

        # Writing package versions to versions.yml
        x = sessionInfo()
        versions <- list(
            R = paste0(x\$R.version\$major, ".", x\$R.version\$minor)
        )
        for(i in seq_along(x\$otherPkgs)){
            pkg <- x\$otherPkgs[[i]]
            versions[[pkg$Package]] <- pkg\$Version
        }
        versions <- list(ORFDEX = versions)
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

}
