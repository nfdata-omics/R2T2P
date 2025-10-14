process R_MERGE_DENOVO_WITH_REF {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/nfdata/riboseqc:v1.2.0-patched"

    input:
    tuple val(meta), path(comp_denovo_gtf)
    path "R-user-lib/*"
    path ref_gtf_Rannotation

    output:
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export R_LIBS_USER=\$PWD/R-user-lib

    Rscript --vanilla ${projectDir}/bin/merge_gtfs.R ${ref_gtf_Rannotation} ${comp_denovo_gtf} FALSE

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
