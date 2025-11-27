process MERGE_TO_ORFQUANT {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/nfdata/riboseqc:v1.3.0-patched"

    input:
    tuple val(meta), path(bam_for_orfquant)

    output:
    tuple val(meta), path("*.RData"), emit: combined_for_orfquant
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def r_files = bam_for_orfquant.collect { file -> "\"${file}\"" }.join(", ")
    """
    Rscript --vanilla -e '
        library(RiboseQC)

        # Running RiboseQC analysis
        RiboseQC:::pool_forORFquant(
            c(${r_files}),
            "${meta.id}.RData"
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
        versions <- list(MERGE_TO_ORFQUANT = versions)
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

    stub:
    """
    touch ${meta.id}.RData

    Rscript -e '
        library(RiboseQC)

        # Writing package versions to versions.yml
        x = sessionInfo()
        versions <- list(
            R = paste0(x\$R.version\$major, ".", x\$R.version\$minor)
        )
        for(i in seq_along(x\$otherPkgs)){
            pkg <- x\$otherPkgs[[i]]
            versions[[pkg$Package]] <- pkg\$Version
        }
        versions <- list(MERGE_TO_ORFQUANT = versions)
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

}
