process ORFQUANT {
    tag "${meta.id}"
    label 'process_high','process_long','process_high_memory'

    container "docker.io/nfdata/orfquant:v481ec99847e9d253a11333a5ae12f0a760338501"

    input:
    tuple val(meta), path(combined_for_orfquant)
    path "R-user-lib/*"
    path gtf_Rannot

    output:
    tuple val(meta), path("*_Detected_ORFs.gtf"),       emit: gtf
    tuple val(meta), path("*_Protein_sequences.fasta"), emit: fasta
    tuple val(meta), path("*_final_ORFquant_results"),  emit: final_results
    tuple val(meta), path("*_tmp_ORFquant_results"),    emit: tmp_results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export R_LIBS_USER=\$PWD/R-user-lib

    Rscript --vanilla -e '
        library("ORFquant")

        run_ORFquant(
            for_ORFquant_file = "${combined_for_orfquant}",
            annotation_file = "${gtf_Rannot}",
            n_cores = ${task.cpus},
            $args
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
        versions <- list(ORFQUANT = versions)
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

    stub:
    """
    Rscript -e '
        library("ORFquant")

        # Writing package versions to versions.yml
        x = sessionInfo()
        versions <- list(
            R = paste0(x\$R.version\$major, ".", x\$R.version\$minor)
        )
        for(i in seq_along(x\$otherPkgs)){
            pkg <- x\$otherPkgs[[i]]
            versions[[pkg$Package]] <- pkg\$Version
        }
        versions <- list(ORFQUANT = versions)
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """
}
