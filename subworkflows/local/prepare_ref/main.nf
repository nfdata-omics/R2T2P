//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF              } from '../../../modules/nf-core/gunzip'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../../modules/nf-core/untar'

include { CUSTOM_GETCHROMSIZES              } from '../../../modules/nf-core/custom/getchromsizes'
include { GFFREAD                           } from '../../../modules/nf-core/gffread'
include { STAR_GENOMEGENERATE               } from '../../../modules/nf-core/star/genomegenerate'
include { FASTATOTWOBIT                     } from '../../../modules/local/fatotwobit'
include { PREPARE_ANNOTATION_FILES          } from '../../../modules/local/prepare_annotation_files'

workflow PREPARE_REF {

    take:
    fasta                    // file: /path/to/genome.fasta (optional!)
    gtf                      // file: /path/to/genome.gtf
    gff                      // file: /path/to/genome.gff
    star_index               // directory: /path/to/star/index/ (optional!)

    main:

    // Versions collector
    ch_versions = Channel.empty()

    // Uncompress GTF or GFF, and convert GFF to GTF
    ch_gtf = Channel.empty()
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            GUNZIP_GTF( gtf.map { [ [:], it ] } )
            ch_gtf      = GUNZIP_GTF.out.gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.value(file(gtf, checkIfExists: true))
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            GUNZIP_GFF( gff.map { [ [:], it ] } )
            ch_gff      = GUNZIP_GFF.out.gunzip
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = gff.map { [ [:], it ] }
        }
        ch_gtf      = GFFREAD(ch_gff, []).gtf.map { it[1] }
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    // Uncompress FASTA if needed
    ch_fasta = Channel.of([])
    if (params.fasta.endsWith('.gz')) {
        GUNZIP_FASTA( fasta.map { [ [:], it ] } )
        ch_fasta    = GUNZIP_FASTA.out.gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = fasta
    }

    // get FASTA index and sizes of the chromosomes
    ch_fai         = Channel.empty()
    ch_chrom_sizes = Channel.empty()

    CUSTOM_GETCHROMSIZES( ch_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    // Build a channel with the STAR index
    ch_star_index = Channel.empty()
    if (params.star_index) {    // If a STAR index is provided, use it
        if (params.star_index.endsWith('.tar.gz')) {
            UNTAR_STAR_INDEX( star_index.map { [ [:], it ] } )
            ch_star_index = UNTAR_STAR_INDEX.out.untar.map { it[1] }
            ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
        } else {
            ch_star_index = star_index
        }
    } else  {         // Otherwise build new STAR index
        STAR_GENOMEGENERATE(
            ch_fasta.map { [ [:], it ] }
        )
        ch_star_index = STAR_GENOMEGENERATE.out.index.map { it[1] }
        ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    // Build a 2bit version of the genome fasta file
    FASTATOTWOBIT( ch_fasta.map { [ [id: it.name], it ] } )
    ch_versions = ch_versions.mix(FASTATOTWOBIT.out.versions)

    // Create annotation files with RiboseQC prepare_annotation_files function
    PREPARE_ANNOTATION_FILES(
        FASTATOTWOBIT.out.genome2bit,
        ch_gtf.map { [ [id: it.name], it ] }
    )
    ch_gtf_Rannot = PREPARE_ANNOTATION_FILES.out.gtf_Rannot.map { it[1] }
    ch_bsgenome = PREPARE_ANNOTATION_FILES.out.bsgenome.map { it[1] }

    emit:
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    fai              = ch_fai                    // channel: path(genome.fai)
    chrom_sizes      = ch_chrom_sizes            // channel: path(genome.sizes)
    star_index       = ch_star_index             // channel: path(star/index/)
    gtf_Rannot       = ch_gtf_Rannot             // channel: path(genome.gtf_Rannot)
    bsgenome         = ch_bsgenome               // channel: path(BSgenome.species.assembly/)
    versions         = ch_versions               // channel: [ versions.yml ]

}
