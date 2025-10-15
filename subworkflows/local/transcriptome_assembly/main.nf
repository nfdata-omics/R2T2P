//
// De novo transcriptome assembly
//

include { SAMTOOLS_MERGE                   } from '../../../modules/nf-core/samtools/merge'
include { STRINGTIE_STRINGTIE as STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie'
include { GAWK as FILTER_UNDEFINED_STRAND  } from '../../../modules/nf-core/gawk'
include { GFFCOMPARE                       } from '../../../modules/nf-core/gffcompare'

include { MERGE_WITH_REF_ANNOTATION         } from '../merge_with_ref_annotation'
include { MERGE_WITH_USER_PROVIDED          } from '../merge_with_user_provided'

workflow TRANSCRIPTOME_ASSEMBLY {
    take:
    ch_bam   // channel: [ val(meta), [ bam ] ]
    ch_fasta // channel: path(fasta)
    ch_fai   // channel: path(fai)
    ch_gtf   // channel: path(gtf)
    ch_bsgenome_dir // channel: path(bsgenome_dir)
    ch_gtf_Rannotation // channel: path(gtf_Rannotation)
    ch_user_gtf // channel: path(gtf) (optional, user-provided annotation)

    main:

    ch_versions = Channel.empty()

    // Merge all those BAM files
    SAMTOOLS_MERGE (
        ch_bam.collect { it[1] }.map { [ ["id": "merged_bams"], it[0] ] }, // get list of bam files
        ch_fasta.map { [ [:], it ] },
        ch_fai.map { [ [:], it ] },
        [[], file("$projectDir/assets/NO_FILE")]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    SAMTOOLS_MERGE.out.bam
        .map { [ it[0] + ["strandedness": "forward"], it[1] ] }
        .set { ch_merged_bam }

    // Run StringTie to assemble transcripts from merged alignments
    STRINGTIE (
        ch_merged_bam,
        ch_gtf
    )

    // Keep only transcripts with strand defined (column 7 ≠ “.”) from StringTie GTF
    FILTER_UNDEFINED_STRAND (
        STRINGTIE.out.transcript_gtf,
        [],
        false
    )

    if ( params.user_provided_annotation ) {
        // If the user provided an annotation GTF, merge StringTie GTF with it
        MERGE_WITH_USER_PROVIDED (
            FILTER_UNDEFINED_STRAND.out.output,
            ch_user_gtf
        )
        ch_versions = ch_versions.mix(MERGE_WITH_USER_PROVIDED.out.versions)
        ch_gtf_to_use = MERGE_WITH_USER_PROVIDED.out.gtf
    } else {
        ch_gtf_to_use = FILTER_UNDEFINED_STRAND.out.output
    }

    // MERGE_WITH_REF_ANNOTATION (
    //     ch_gtf_to_use,
    //     ch_gtf,
    //     ch_bsgenome_dir,
    //     ch_gtf_Rannotation
    // )

    emit:
    // bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
