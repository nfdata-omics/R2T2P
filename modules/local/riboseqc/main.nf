process RIBOSEQC {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/nfdata/riboseqc:v1.3.0-patched"

    input:
    tuple val(meta), path(bam)
    path "R-user-lib/*"
    path gtf_Rannot

    output:
    tuple val(meta), path("*.RData"), emit: riboseqc_rdata
    tuple val(meta), path("*.bedgraph"), emit: bedgraph
    tuple val(meta), path("*_strandedness"), emit: strandedness
    tuple val(meta), path("*_counts_regions"), emit: counts_regions
    tuple val(meta), path("*.bam_for_ORFquant"), emit: bam_for_orfquant, optional: true
    tuple val(meta), path("*.bam_results_RiboseQC"), emit: bam_results_riboseqc, optional: true
    tuple val(meta), path("*.bam_results_RiboseQC_all"), emit: bam_results_riboseqc_all, optional: true
    tuple val(meta), path("*.bam_P_sites_calcs"), emit: bam_p_sites_calcs, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export R_LIBS_USER=\$PWD/R-user-lib

    Rscript --vanilla -e '
        library(RiboseQC)
        library(msa)

        # Loading RiboseQC annotation
        load_annotation("${gtf_Rannot}")

        # Running RiboseQC analysis
        RiboseQC_analysis(
            annotation_file = "${gtf_Rannot}",
            bam_file = "${bam}",
            ${args}
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
        versions <- list(RIBOSEQC = versions)
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

    stub:
    """
    touch ${meta.id}_strandedness
    touch ${meta.id}_counts_regions
    touch ${meta.id}.RData
    touch ${meta.id}.bedgraph
    touch ${meta.id}.bam_for_ORFquant
    touch ${meta.id}.bam_results_RiboseQC
    touch ${meta.id}.bam_results_RiboseQC_all
    touch ${meta.id}.bam_P_sites_calcs

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
        versions <- list(RIBOSEQC = versions)
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

}
