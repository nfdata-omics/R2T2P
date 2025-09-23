//
// De novo transcriptome assembly
//

include { GAWK as FILTER_UNDEFINED_STRAND  } from '../../../modules/nf-core/gawk'
include { GFFCOMPARE                       } from '../../../modules/nf-core/gffcompare'
include { R_MERGE_DENOVO_WITH_REF          } from '../../../modules/local/r_merge_denovo_with_ref'

workflow MERGE_WITH_REF_ANNOTATION {
    take:
    ch_merged_gtf  // channel: [ val(meta), [ gtf ] ]
    ch_ref_gtf     // channel: gtf
    ch_bsgenome_dir // channel: path(bsgenome_dir)
    ch_gtf_Rannotation // channel: path(gtf_Rannotation)

    main:

    ch_versions = Channel.empty()

    // Compare (denovo vs ) the merged GTF with user annotation
    GFFCOMPARE (
        ch_merged_gtf.map { [ ["id": "comp_denovo"], it[1] ] },
        [[], [], []], // no fasta
        ch_ref_gtf.map { [ [:], it ] }
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

//   Rscript $1"/e_merge_gtfs.R" $6 comp_denovo.annotated_ok.gtf FALSE $1

    // Concatenate the newly annotated transcripts to the reference GTF
//   cat comp_denovo.annotated_ok.gtf $4 > combined_new_annotated.gtf

    emit:
    // bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
