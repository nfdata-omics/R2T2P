process STAR_GENOMEGENERATE {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "docker.io/nfdata/star:v2.7.9a"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("star")  , emit: index
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def args_list   = args.tokenize()
    def memory      = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    if (args_list.contains('--genomeSAindexNbases')) {
        """
        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta \\
            --runThreadN $task.cpus \\
            $memory \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """
    } else {
        """
        samtools faidx $fasta
        NUM_BASES=`gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${fasta}.fai`

        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta \\
            --runThreadN $task.cpus \\
            --genomeSAindexNbases \$NUM_BASES \\
            $memory \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """
    }

    stub:
        """
        mkdir star
        touch star/Genome
        touch star/Log.out
        touch star/SA
        touch star/SAindex
        touch star/chrLength.txt
        touch star/chrName.txt
        touch star/chrNameLength.txt
        touch star/chrStart.txt
        touch star/genomeParameters.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """
}
