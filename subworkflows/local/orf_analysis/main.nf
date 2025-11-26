//
// ORF quantification and differential analysis
//

include { MERGE_TO_ORFQUANT  } from '../../../modules/local/merge_to_orfquant/main'
include { ORFQUANT           } from '../../../modules/local/orfquant/main'
include { ORFDEX             } from '../../../modules/local/orfdex/main'

workflow ORF_ANALYSIS {

    take:
    ch_bam                   // channel: [ val(meta), path(bam) ]
    ch_bam_for_orfquant      // channel: [ val(meta), path(bam_for_orfquant) ]
    ch_bsgenome              // channel: bsgenome
    ch_gtf_Rannot            // channel: gtf_Rannot

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
        ch_bsgenome,
        ch_gtf_Rannot
    )
    ch_versions = ch_versions.mix(ORFQUANT.out.versions)

    if ( params.control_label ) {
        // Filter to conditions that are not the control, extract unique condition names,
        // and create contrast metadata
        ch_bam
            .filter { meta, _file -> meta.condition != params.control_label }
            .map { meta, _file -> [ meta.condition ]}
            .unique()
            .map { contrast -> [id: "${contrast.join(",")}_vs_${params.control_label}",
                condition: contrast, control: [params.control_label] ] }
            .set { ch_conditions }
    } else {
        // If no control label is provided, set ch_conditions to empty channel to skip DE analysis
        ch_conditions = channel.empty()
    }

    // Combine contrasts with all count files and filter to only relevant files for each contrast
    ch_conditions
        .combine( ch_bam )
        .filter { meta_contrast, meta_file, _file ->
            meta_file.condition in meta_contrast.control + meta_contrast.condition }
        .set { ch_contrasts_full_table }

    // Collect files into contrast-specific tables with metadata columns, ordered by library type
    ch_contrasts_full_table
        .groupTuple(by: 0) // Group by meta_contrast
        .map { meta_contrast, meta_files, files ->
            // Combine and sort the data
            [meta_files, files].transpose()
                .sort { a, b ->
                    def order = [ "RNA": 0, "Ribo": 1 ]
                    def cmp = order.get(a[0].library_type, 99) <=> order.get(b[0].library_type, 99)
                    if (cmp != 0) return cmp
                    return a[0].id <=> b[0].id
                }
                .collect { meta_file, file -> [meta_contrast, meta_file, file] }
        }
    .flatten()
    .collate(3) // Group back into tuples of 3 elements
    .collectFile( { meta_contrast, meta_file, file ->
       [ "${meta_contrast.id}.txt", [ file.name, meta_file.library_type, meta_file.condition,
            (meta_file.condition in meta_contrast.control ? "TRUE" : "FALSE" ) ].join('\t') + "\n" ]
        }, sort: "index")
    .map { file -> [ file.name.replace('.txt', ''), file ]  }
    .set { ch_contrasts_table }

    // Group files by contrast, merge with table metadata, and prepare final channel
    ch_contrasts_full_table
        .map {  meta_contrast, _meta_file, file -> [ meta_contrast, file ]}
        .groupTuple()
        .map { meta, files -> [ meta.id, meta, files ] }
        .join( ch_contrasts_table )
        .map { _id, meta, files, table -> [ meta, files, table ] }
        .set { ch_contrast_table_with_files }

    //
    // Differential ORF analysis
    //
    ORFDEX (
        ch_contrast_table_with_files,
        ORFQUANT.out.final_results.map { _meta, file -> file }
    )
    ch_versions = ch_versions.mix(ORFDEX.out.versions.first())

    emit:
    versions   = ch_versions                              // channel: [ versions.yml ]

}
