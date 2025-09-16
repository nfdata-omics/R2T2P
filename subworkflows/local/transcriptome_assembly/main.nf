
//
// De novo transcriptome assembly
//

include { SAMTOOLS_MERGE } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_VIEW } from '../../../modules/nf-core/samtools/view/main'
include { STRINGTIE_STRINGTIE as STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'

workflow TRANSCRIPTOME_ASSEMBLY {
    take:
    ch_bam   // channel: [ val(meta), [ bam ] ]
    ch_fasta // channel: [ val(meta), fasta ]
    ch_fai   // channel: [ val(meta), fai ]
    ch_gtf   // channel: [ val(meta), gtf ]

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_MERGE (
        ch_bam,
        ch_fasta,
        ch_fai,
        [[], file("$projectDir/assets/NO_FILE")]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    // SAMTOOLS_VIEW (
    //     SAMTOOLS_MERGE.out.bam.join(file("")),
    //     ch_fasta,
    //     file(""),
    //     'bai'
    // )

    SAMTOOLS_MERGE.out.bam
        .map { [ it[0] + ["strandedness": "forward"], it[1] ] }
        .set { ch_merged_bam }

    STRINGTIE (
        ch_merged_bam,
        ch_gtf
    )

    emit:
    // bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
