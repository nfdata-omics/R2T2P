//
// Differential analysis of multiple conditions vs a baseline
//

include { MULTI_DE } from '../../../modules/local/multi_de/main'
include { GET_GO_PACKAGE } from '../../../modules/local/get_go_package/main'

workflow DIFFERENTIAL_ANALYSIS {
    take:
    ch_counts_regions        // channel: [ val(meta), path(counts_regions) ]
    ch_bsgenome              // channel: bsgenome
    ch_gtf_Rannot            // channel: gtf_Rannot

    main:
    ch_versions = channel.empty()

    if ( params.control_label ) {
        // Filter to conditions that are not the control, extract unique condition names,
        // and create contrast metadata
        ch_counts_regions
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
        .combine( ch_counts_regions )
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
                    def cmp = order.get(a[0].assay_type, 99) <=> order.get(b[0].assay_type, 99)
                    if (cmp != 0) return cmp
                    return a[0].id <=> b[0].id
                }
                .collect { meta_file, file -> [meta_contrast, meta_file, file] }
        }
    .flatten()
    .collate(3) // Group back into tuples of 3 elements
    .collectFile( { meta_contrast, meta_file, file ->
       [ "${meta_contrast.id}.txt", [ file.name, meta_file.assay_type, meta_file.condition,
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
    // Get the org.db package name if GO analysis is requested
    //
    if ( params.go_package ) {
        GET_GO_PACKAGE ( params.go_package )
        org_db_package = GET_GO_PACKAGE.out.org_db_package.collect()
    } else {
        org_db_package = file("$projectDir/assets/NO_FILE")
    }

    //
    // Run multi-DE analysis for each contrast
    //
    MULTI_DE (
        ch_contrast_table_with_files,
        ch_bsgenome,
        ch_gtf_Rannot,
        params.go_package,
        org_db_package
    )
    ch_versions = ch_versions.mix(MULTI_DE.out.versions.first())

    emit:
    versions  = ch_versions                         // channel: [ versions.yml ]

}
