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

    ch_versions = channel.empty()
    ch_gff_stats = channel.empty()

    // Collect BAM files and create merged meta validating that strandedness is consistent across all BAM files
    ch_bam
        .map { meta, bam -> [meta.strandedness, bam] }
        .collect(flat: false)
        .map { collected_items ->
            // Extract strandedness values and BAM files
            def unique_strandedness = collected_items.collect { item -> item[0] }.unique()
            def bam_files = collected_items.collect { item -> item[1] }
            // Check if all strandedness values are the same
            if (unique_strandedness.size() > 1) {
                error "Inconsistent strandedness values found: ${unique_strandedness}. All BAM files must have the same strandedness."
            }
            // Create new meta with merged information
            def merged_meta = [ "id": "merged_bams", "strandedness": unique_strandedness[0] ]
            return [merged_meta, bam_files]
        }
        .set { ch_collected_bams }

    // Merge all those BAM files
    SAMTOOLS_MERGE (
        ch_collected_bams,
        ch_fasta.map { file -> [ [:], file ] },
        ch_fai.map { file -> [ [:], file ] },
        [[], []]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    // Run StringTie to assemble transcripts from merged alignments
    STRINGTIE (
        SAMTOOLS_MERGE.out.bam,
        ch_gtf
    )
    ch_versions = ch_versions.mix(STRINGTIE.out.versions)

    // Keep only transcripts with strand defined (column 7 ≠ “.”) from StringTie GTF
    FILTER_UNDEFINED_STRAND (
        STRINGTIE.out.transcript_gtf,
        [],
        false
    )
    ch_versions = ch_versions.mix(FILTER_UNDEFINED_STRAND.out.versions)

    if ( params.user_provided_annotation ) {
        // If the user provided an annotation GTF, merge StringTie GTF with it
        MERGE_WITH_USER_PROVIDED (
            FILTER_UNDEFINED_STRAND.out.output,
            ch_user_gtf
        )
        ch_versions = ch_versions.mix(MERGE_WITH_USER_PROVIDED.out.versions)
        ch_gtf_to_use = MERGE_WITH_USER_PROVIDED.out.gtf
        ch_gff_stats = ch_gff_stats.mix(MERGE_WITH_USER_PROVIDED.out.gffcompare_stats)
    } else {
        ch_gtf_to_use = FILTER_UNDEFINED_STRAND.out.output.map { _meta, file -> [ file ] }
    }

    MERGE_WITH_REF_ANNOTATION (
        ch_gtf_to_use,
        ch_gtf,
        ch_bsgenome_dir,
        ch_gtf_Rannotation
    )
    ch_gff_stats = ch_gff_stats.mix(MERGE_WITH_REF_ANNOTATION.out.gffcompare_stats)
    ch_versions = ch_versions.mix(MERGE_WITH_REF_ANNOTATION.out.versions)

    emit:
    gtf        = MERGE_WITH_REF_ANNOTATION.out.gtf        // channel: [ val(meta), [ gtf ] ]
    gtf_Rannot = MERGE_WITH_REF_ANNOTATION.out.gtf_Rannot // channel: [ val(meta), [ gtf_Rannot ] ]
    gff_stats  = ch_gff_stats                             // channel: [ stats ]
    versions   = ch_versions                              // channel: [ versions.yml ]
}
