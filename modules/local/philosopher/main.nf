process PHILOSOPHER {
    tag "${meta.id}"
    label 'process_low'

    container "docker.io/fcyucn/fragpipe:23.1"

    input:
    tuple val(meta), path(db_fasta)

    output:
    tuple val(meta), path("*.fasta.fas"), emit: fasta_with_decoys_and_contam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    PHILO_EXE="/fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/Philosopher/philosopher-v5.1.2"
    \$PHILO_EXE workspace --init
    \$PHILO_EXE database --custom ${db_fasta} --contam
    \$PHILO_EXE workspace --clean

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        philosopher: \$(\$PHILO_EXE version 2>&1 | sed -n "s/.*version=\\(v[0-9.]*\\).*/\\1/p")
    END_VERSIONS
    """

    stub:
    """
    touch \$(date "+%F")-decoys-contam-${meta.id}.fasta.fas

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        philosopher: \$(\$PHILO_EXE version 2>&1 | sed -n "s/.*version=\\(v[0-9.]*\\).*/\\1/p")
    END_VERSIONS
    """

}
