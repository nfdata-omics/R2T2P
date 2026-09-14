//
// Create and analyze proteomics database
//

include { CREATE_PROTEIN_DB_WRITE_DB      } from '../../../modules/local/create_protein_db_write_db/main'
include { PHILOSOPHER                     } from '../../../modules/local/philosopher/main'
include { CREATE_PROTEIN_DB_REPLACE_NAMES } from '../../../modules/local/create_protein_db_replace_names/main'
include { FRAGPIPE                        } from '../../../modules/local/fragpipe/main'

workflow PROTEOMICS {

    take:
    ch_bsgenome              // channel: bsgenome
    ch_gtf_Rannot            // channel: gtf_Rannot (ref or product of R module, depending)
    ch_orfquant_fasta        // channel: fasta with proteins sequences from ORFquant
    workflow_file            // value channel: path to fragpipe workflow file
    ch_tools_folder          // value channel: folder with the external tools for fragpipe
    ch_diann_folder          // value channel: folder with the diann installation for fragpipe
    ch_manifest              // value channel: path to fragpipe manifest file
    ch_annotation_file       // value channel: path to tmt_annotation file (optional)

    main:

    ch_versions = channel.empty()

    // skip proteomics workflow if no fragpipe manifest is provided
    if ( params.fragpipe_manifest == null ) {
        ch_protein_dbs = channel.empty()
    } else {
        ch_protein_dbs = ch_orfquant_fasta
    }

    CREATE_PROTEIN_DB_WRITE_DB (
        ch_bsgenome,
        ch_gtf_Rannot,
        ch_protein_dbs
    )
    ch_versions = ch_versions.mix(CREATE_PROTEIN_DB_WRITE_DB.out.versions)

    CREATE_PROTEIN_DB_WRITE_DB.out.annot_prot_fasta.map{ file -> [["id": file.baseName], file] }
        .concat( CREATE_PROTEIN_DB_WRITE_DB.out.orfquant_prot_fasta.map{ file -> [["id": file.baseName], file] }  )
        .concat( CREATE_PROTEIN_DB_WRITE_DB.out.annot_and_orfquant_prot_fasta.map{ file -> [["id": file.baseName], file] }  )
        .set{ ch_protein_dbs }

    PHILOSOPHER (
        ch_protein_dbs
    )
    ch_versions = ch_versions.mix(PHILOSOPHER.out.versions.first())

    ch_protein_dbs
        .join( PHILOSOPHER.out.fasta_with_decoys_and_contam )
        .set { ch_protein_db_and_fasta_with_decoys_and_contam }

    CREATE_PROTEIN_DB_REPLACE_NAMES (
        ch_protein_db_and_fasta_with_decoys_and_contam
    )
    ch_versions = ch_versions.mix(CREATE_PROTEIN_DB_REPLACE_NAMES.out.versions.first())

    ch_manifest
        .splitCsv( header: ["path", "experiment_name", "bioreplicate", "data_type"], sep: '\t' )
        .map { row -> file(row.path) }
        .collect()
        .set { ch_fragpipe_mzml_files }

    FRAGPIPE (
        CREATE_PROTEIN_DB_REPLACE_NAMES.out.fasta_with_decoys_and_contam,
        workflow_file,
        ch_manifest,
        ch_fragpipe_mzml_files,
        ch_annotation_file,
        ch_tools_folder,
        ch_diann_folder
    )
    ch_versions = ch_versions.mix(FRAGPIPE.out.versions.first())

    emit:
    versions   = ch_versions                              // channel: [ versions.yml ]

}
