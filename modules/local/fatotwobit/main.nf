process FASTATOTWOBIT {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/nfdata/fatotwobit:v2025-09-10"

    input:
    tuple val(meta), path(genome_fasta)

    output:
    tuple val(meta), path("*.2bit"), emit: genome2bit
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    faToTwoBit ${genome_fasta} ${prefix}.2bit

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FaToTwoBit: v2025-09-10
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.2bit

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FaToTwoBit: v2025-09-10
    END_VERSIONS
    """
}
