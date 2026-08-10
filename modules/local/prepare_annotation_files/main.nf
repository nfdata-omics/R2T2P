process PREPARE_ANNOTATION_FILES {
    tag "${meta2.id}"
    label 'process_single'

    container "docker.io/nfdata/riboseqc:v1.3.0-patched"

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
         library(GenomeInfoDb)
         library(RiboseQC)
         
         find_grch38_ncbi_dir <- function(
             assembly_accession,
             assembly_name = NA_character_
         ) {
             if (!identical(assembly_accession, "GCF_000001405.26")) {
                 stop(
                     "Temporary workaround supports only GCF_000001405.26; received: ",
                     assembly_accession
                 )
             }
         
             prefix <- "GCF_000001405.26_GRCh38"
         
             directory <- paste0(
                 "https://ftp.ncbi.nlm.nih.gov/genomes/all/",
                 "GCF/000/001/405/",
                 prefix
             )
         
             c(directory, prefix)
         }
         
         assignInNamespace(
             "find_NCBI_assembly_ftp_dir",
             find_grch38_ncbi_dir,
             ns = "GenomeInfoDb"
         )
         
         stopifnot(
             nrow(
                 GenomeInfoDb::fetch_assembly_report(
                     "GCF_000001405.26",
                     assembly_name = "GRCh38"
                 )
             ) > 0L
        )
        
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
        versions <- list(PREPARE_ANNOTATION_FILES = versions)
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
        versions <- list(PREPARE_ANNOTATION_FILES = versions)
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

}
