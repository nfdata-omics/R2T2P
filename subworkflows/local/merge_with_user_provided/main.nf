//
// De novo transcriptome assembly
//

include { GAWK as FILTER_UNDEFINED_STRAND  } from '../../../modules/nf-core/gawk'
include { GFFCOMPARE                       } from '../../../modules/nf-core/gffcompare'
include { R_MERGE_DENOVO_WITH_USER         } from '../../../modules/local/r_merge_denovo_with_user'

workflow MERGE_WITH_USER_PROVIDED {
    take:
    ch_merged_gtf  // channel: [ val(meta), [ gtf ] ]
    ch_user_gtf     // channel: gtf

    main:

    ch_versions = channel.empty()

    // Compare the user provided annotation with de novo−only GTF
    GFFCOMPARE (
        ch_user_gtf.map { gtf -> [ ["id": "comp_strg_user"], gtf ] },
        [[], [], []], // no fasta
        ch_merged_gtf.map { _meta, gtf -> [ [:], gtf ] }
    )
    ch_versions = ch_versions.mix(GFFCOMPARE.out.versions)

    // Keep only transcripts with strand defined (column 7 ≠ “.”)
    FILTER_UNDEFINED_STRAND (
        GFFCOMPARE.out.annotated_gtf,
        [],
        false
    )
    ch_versions = ch_versions.mix(FILTER_UNDEFINED_STRAND.out.versions)

    // Merge StringTie GTF and user annotation GTF via an R script
    R_MERGE_DENOVO_WITH_USER (
        FILTER_UNDEFINED_STRAND.out.output,
        ch_merged_gtf.map { _meta, gtf -> [ gtf ] },
        ch_user_gtf
    )
    ch_versions = ch_versions.mix(R_MERGE_DENOVO_WITH_USER.out.versions)

    emit:
    gtf              = R_MERGE_DENOVO_WITH_USER.out.gtf        // channel: [ gft ]
    gffcompare_stats = GFFCOMPARE.out.stats.map { _meta, stats -> [ stats ] }  // channel: [ stats ]
    versions         = ch_versions                             // channel: [ versions.yml ]
}
