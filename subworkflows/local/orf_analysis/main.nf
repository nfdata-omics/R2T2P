//
// ORF quantification and differential analysis
//

include { MERGE_TO_ORFQUANT  } from '../../../modules/local/merge_to_orfquant/main'
include { ORFQUANT           } from '../../../modules/local/orfquant/main'

workflow ORF_ANALYSIS {

    take:
    ch_bam_for_orfquant // channel: [ val(meta), path(bam_for_orfquant) ]
    ch_gtf_Rannot       // channel: gtf_Rannot

    main:

    ch_versions = channel.empty()

    //
    // Merge files for ORFquant
    //

    MERGE_TO_ORFQUANT (
        ch_bam_for_orfquant
            .collect { _meta, file -> file }
            .map { files -> [ ["id": "combined_forORFquant"], files ] }
    )
    ch_versions = ch_versions.mix(MERGE_TO_ORFQUANT.out.versions)

    //
    // ORF quantification with ORFquant
    //

    ORFQUANT (
        MERGE_TO_ORFQUANT.out.combined_for_orfquant,
        ch_gtf_Rannot
    )
    ch_versions = ch_versions.mix(ORFQUANT.out.versions)

    //
    // Differential ORF analysis
    //


    emit:
    versions   = ch_versions                              // channel: [ versions.yml ]

}
