process DIFF_PROTEOMICS {
    tag "${meta.id}"
    label 'process_single', 'process_high_memory'

    container "docker.io/nfdata/riboseqc:v2.1.0-patched"

    input:
    tuple val(meta), path(search_folder), path(tmt_annotation_files, stageAs: 'tmt_annotations??/*')
    path manifest_file
    path "R-user-lib/*"
    path gtf_Rannot
    path orfquant_results
    val data_type
    val control_label

    output:
    path "*.RData", emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def annotation_files = data_type == "DDA_TMT" ? "--annotation-files \"" + [tmt_annotation_files].flatten().join(',') + "\"" : ''
    """
    export R_LIBS_USER=\$PWD/R-user-lib

    Rscript --vanilla \"${projectDir}/bin/diff_proteomics.R\" ${args} \
        --data-type \"${data_type}\" \
        --gtf-rannot \"${gtf_Rannot}\" \
        --orfquant-results \"${orfquant_results}\" \
        --search-fold \"${search_folder}\" \
        --manifest \"${manifest_file}\" \
        --baseline \"${control_label}\" \
        ${annotation_files}

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
        versions <- list(DIFF_PROTEOMICS = versions)
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """


    stub:
    """

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
        versions <- list(DIFF_PROTEOMICS = versions)
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

}
