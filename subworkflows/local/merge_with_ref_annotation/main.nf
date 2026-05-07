//
// De novo transcriptome assembly
//

include { GAWK as FILTER_UNDEFINED_STRAND  } from '../../../modules/nf-core/gawk'
include { GFFCOMPARE                       } from '../../../modules/nf-core/gffcompare'
include { R_MERGE_DENOVO_WITH_REF          } from '../../../modules/local/r_merge_denovo_with_ref'
include { CAT_GTF                          } from '../../../modules/local/cat_gtf'

workflow MERGE_WITH_REF_ANNOTATION {
    take:
    ch_merged_gtf  // channel: gtf
    ch_ref_gtf     // channel: gtf
    ch_bsgenome_dir // channel: path(bsgenome_dir)
    ch_gtf_Rannotation // channel: path(gtf_Rannotation)

    main:

    ch_versions = channel.empty()

    // Compare (denovo vs ) the merged GTF with user annotation
    GFFCOMPARE (
        ch_merged_gtf.map { gtf -> [ ["id": "comp_denovo"], gtf ] },
        [[], [], []], // no fasta
        ch_ref_gtf.map { gtf -> [ [:], gtf ] }
    )
    ch_versions = ch_versions.mix(GFFCOMPARE.out.versions)

    // Keep only transcripts with strand defined (column 7 ≠ “.”)
    FILTER_UNDEFINED_STRAND (
        GFFCOMPARE.out.annotated_gtf,
        [],
        false
    )
    ch_versions = ch_versions.mix(FILTER_UNDEFINED_STRAND.out.versions)

    // Combine denovo with reference GTF via an R script
    R_MERGE_DENOVO_WITH_REF (
        FILTER_UNDEFINED_STRAND.out.output,
        ch_bsgenome_dir,
        ch_gtf_Rannotation
    )
    ch_versions = ch_versions.mix(R_MERGE_DENOVO_WITH_REF.out.versions)

    // Concatenate the newly annotated transcripts to the reference GTF
    CAT_GTF (
        R_MERGE_DENOVO_WITH_REF.out.gtf,
        ch_ref_gtf
    )
    ch_versions = ch_versions.mix(CAT_GTF.out.versions)

    emit:
    gtf              = CAT_GTF.out.gtf                         // channel: [ val(meta), [ gtf ] ]
    gtf_Rannot       = R_MERGE_DENOVO_WITH_REF.out.gtf_Rannot  // channel: [ gtf_Rannot ]
    gffcompare_stats = GFFCOMPARE.out.stats.map { _meta, stats -> [ stats ] }  // channel: [ stats ]
    versions         = ch_versions                             // channel: [ versions.yml ]
}
