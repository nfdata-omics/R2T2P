//
// Two-pass alignment with STAR
//

include { STAR_ALIGN as STAR_FIRST_ALIGN                  } from '../../../modules/nf-core/star/align/main'
include { SAMTOOLS_INDEX as FIRST_SAMTOOLS_INDEX          } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS as FIRST_BAM_STATS_SAMTOOLS  } from '../../../subworkflows/nf-core/bam_stats_samtools/main'
include { CREATE_FIRSTPASS_JUNCTIONS                      } from '../../../modules/local/create_firstpass_junctions/main'
include { STAR_ALIGN as STAR_WITH_NOVEL_JUNCT             } from '../../../modules/nf-core/star/align/main'
include { SAMTOOLS_INDEX as SECOND_SAMTOOLS_INDEX         } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS as SECOND_BAM_STATS_SAMTOOLS } from '../../../subworkflows/nf-core/bam_stats_samtools/main'
include { RIBOSEQC                                        } from '../../../modules/local/riboseqc/main'
include { UCSC_BEDGRAPHTOBIGWIG                           } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow TWO_PASS_ALIGNMENT {
    take:
    ch_reads        // channel: [ val(meta), [ fastq ] ]
    ch_star_index   // path(genome_index)
    ch_gtf          // path(gtf)
    ch_fasta        // path(fasta)
    ch_fai          // path(fai)
    ch_chrom_sizes  // channel: chrom_sizes
    ch_samplesheet  // samplesheet channel for ordering
    ch_bsgenome     // path(bsgenome)
    ch_gtf_Rannot   // path(gtf_Rannot)

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // Map reads with STAR
    //
    STAR_FIRST_ALIGN (
        ch_reads,
        ch_star_index.map { file -> [ [:], file ] },
        ch_gtf.map { file -> [ [:], file ] },
        "$projectDir/assets/NO_FILE", // empty arguments for additional_junctions
        false, // star_ignore_sjdbgtf
        "", // seq_platform
        "" // seq_center
    )
    ch_versions = ch_versions.mix(STAR_FIRST_ALIGN.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(STAR_FIRST_ALIGN.out.log_final.collect{ _meta, file -> file })

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //

    FIRST_SAMTOOLS_INDEX ( STAR_FIRST_ALIGN.out.bam_sorted_aligned )
    ch_versions = ch_versions.mix(FIRST_SAMTOOLS_INDEX.out.versions.first())

    STAR_FIRST_ALIGN.out.bam_sorted_aligned
        .join(FIRST_SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(FIRST_SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }

    FIRST_BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta.map { file -> [ [:], file ] }.combine(ch_fai).collect() )

    ch_multiqc_files  = ch_multiqc_files.mix( FIRST_BAM_STATS_SAMTOOLS.out.stats.collect{ _meta, file -> file } )
        .mix( FIRST_BAM_STATS_SAMTOOLS.out.flagstat.collect{ _meta, file -> file } )
        .mix( FIRST_BAM_STATS_SAMTOOLS.out.idxstats.collect{ _meta, file -> file } )

    //
    // Merge all the novel junction tables into a single file
    //
    ch_samplesheet
        .toList()
        .flatMap { list ->
            list.withIndex().collect { item, index ->
                [item[0], index]  // [meta, position]
            }
        }
        .join(STAR_FIRST_ALIGN.out.pass1_spl_juc_tab)
        .toSortedList { a, b -> a[1] <=> b[1] }
        .filter { list -> !list.isEmpty() }
        .map { list -> [["id": "merged_junctions"], list.collect { _meta, _position, junction -> junction }] }
        .set { ch_merged_junctions }

    //
    // Get the table of novel junctions from the first STAR alignment
    //
    CREATE_FIRSTPASS_JUNCTIONS(
        ch_merged_junctions,
        ch_bsgenome,
        ch_gtf_Rannot
    )
    ch_versions = ch_versions.mix(CREATE_FIRSTPASS_JUNCTIONS.out.versions)

    // Match the reads wit the corresponding junction table
    ch_reads
        .combine( CREATE_FIRSTPASS_JUNCTIONS.out.pass1_junctions.map { _meta, file -> file } )
        .set { ch_reads_with_junctions }
    // Split into two channels:
    ch_reads_ordered = ch_reads_with_junctions.map { meta, reads, _junctions -> [meta, reads] }
    ch_junctions_ordered = ch_reads_with_junctions.map { _meta, _reads, junctions -> junctions }

    //
    // Second STAR alignment using novel junctions
    //
    STAR_WITH_NOVEL_JUNCT (
        ch_reads_ordered,
        ch_star_index.map { file -> [ [:], file ] },
        ch_gtf.map { file -> [ [:], file ] },
        ch_junctions_ordered, // channel for additional junctions
        false, // star_ignore_sjdbgtf
        "", // seq_platform
        "" // seq_center
    )
    ch_versions = ch_versions.mix(STAR_WITH_NOVEL_JUNCT.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(STAR_WITH_NOVEL_JUNCT.out.log_final.collect{ _meta, file -> file })

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //

    SECOND_SAMTOOLS_INDEX ( STAR_WITH_NOVEL_JUNCT.out.bam_sorted_aligned )
    ch_versions = ch_versions.mix( SECOND_SAMTOOLS_INDEX.out.versions.first() )

    STAR_WITH_NOVEL_JUNCT.out.bam_sorted_aligned
        .join(SECOND_SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SECOND_SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }

    SECOND_BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta.map { file -> [ [:], file ] }.combine(ch_fai).collect() )

    ch_multiqc_files  = ch_multiqc_files.mix( SECOND_BAM_STATS_SAMTOOLS.out.stats.collect{ _meta, file -> file } )
        .mix( SECOND_BAM_STATS_SAMTOOLS.out.flagstat.collect{ _meta, file -> file } )
        .mix( SECOND_BAM_STATS_SAMTOOLS.out.idxstats.collect{ _meta, file -> file } )

    //
    // Quality control with Ribo-seQC
    //

    RIBOSEQC (
        STAR_WITH_NOVEL_JUNCT.out.bam_sorted_aligned,
        ch_bsgenome,
        ch_gtf_Rannot
    )
    ch_versions = ch_versions.mix(RIBOSEQC.out.versions.first())

    //
    // Convert bedGraph to bigWig
    //

    UCSC_BEDGRAPHTOBIGWIG (
        RIBOSEQC.out.bedgraph
            .flatMap { meta, files -> files.collect { file -> [ meta, file ] } },
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())

    emit:
    bam           = STAR_WITH_NOVEL_JUNCT.out.bam_sorted_aligned // channel: [ val(meta), [ bam ] ]
    strandedness  = RIBOSEQC.out.strandedness                    // channel: [ val(meta), path(strandedness) ]
    versions      = ch_versions                                  // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files                             // channel: [ stats ]

}
