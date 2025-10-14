process R_MERGE_DENOVO_WITH_USER {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/nfdata/riboseqc:v1.2.0-patched"

    input:
    tuple val(meta), path(comp_strg_user_gft)
    path merged_gtf
    path user_gft

    output:
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export R_LIBS_USER=\$PWD/R-user-lib

    Rscript --vanilla ${projectDir}/bin/merge_strg_user_gtfs.R ${merged_gtf} ${comp_strg_user_gft} ${user_gft} TRUE

    Rscript --vanilla -e '
        library(RiboseQC)

        # Writing package versions to versions.yml
        x = sessionInfo()
        versions <- list(
            R = paste0(x\$R.version\$major, ".", x\$R.version\$minor)
        )
        for(i in seq_along(x\$otherPkgs)){
            pkg <- x\$otherPkgs[[i]]
            versions[[pkg\$Package]] <- pkg\$Version
        }
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

    stub:
    """
    touch pippo

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
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

}
