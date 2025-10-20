process PREPARE_ANNOTATION_FILES {
    tag "${meta2.id}"
    label 'process_single'

    container "docker.io/nfdata/riboseqc:v1.2.0-patched"

    input:
    tuple val(meta), path(genome_2bit)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta2), path("annotation/*.gtf_Rannot"), emit: gtf_Rannot
    tuple val(meta), path("BSgenome.any.species.GRCh38/"), emit: bsgenome
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export R_LIBS_USER=\$PWD

    Rscript --vanilla -e '
        library(RiboseQC)
        RiboseQC::prepare_annotation_files(
            "annotation",
            "${genome_2bit}",
            "${gtf}",
            "any.species",
            "GRCh38"
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
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

    stub:
    """
    touch annotation/genome.gtf_Rannot
    mkdir -p BSgenome.species.assembly/

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
