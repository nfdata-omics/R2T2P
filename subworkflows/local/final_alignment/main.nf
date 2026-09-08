//
// Final alignment against the new transcriptome
//

include { STAR_ALIGN as STAR_ALIGN_RNA       } from '../../../modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_RIBO      } from '../../../modules/nf-core/star/align/main'
include { SAMTOOLS_INDEX                     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS                 } from '../../../subworkflows/nf-core/bam_stats_samtools/main'
include { RIBOSEQC as RIBOSEQC_RNA           } from '../../../modules/local/riboseqc/main'
include { RIBOSEQC as RIBOSEQC_RIBO          } from '../../../modules/local/riboseqc/main'
include { UCSC_BEDGRAPHTOBIGWIG              } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

    workflow FINAL_ALIGNMENT {
    take:
    ch_reads        // channel: [ val(meta), [ fastq ] ]
    ch_genome_index // channel: genome_index
    ch_gtf          // channel: gtf
    ch_fasta        // channel: fasta
    ch_fai          // channel: fai
    ch_chrom_sizes  // channel: chrom_sizes
    ch_bsgenome     // channel: bsgenome
    ch_gtf_Rannot   // channel: gtf_Rannot

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_reads
        .branch { meta, _fastq ->
            rna:  meta.assay_type == "RNA"
            ribo: meta.assay_type == "Ribo"
        }
    .set { ch_reads_by_type }

    //
    // align RNA reads with STAR
    //

    STAR_ALIGN_RNA (
        ch_reads_by_type.rna,
        ch_genome_index.map { file -> [ [:], file ] },
        ch_gtf.map { file -> [ [:], file ] },
        "$projectDir/assets/NO_FILE", // empty arguments for additional_junctions
        false, // star_ignore_sjdbgtf
        "", // seq_platform
        "" // seq_center
    )
    ch_versions = ch_versions.mix(STAR_ALIGN_RNA.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN_RNA.out.log_final.collect{ _meta, log -> log })

    //
    // Quality control with Ribo-seQC for RNAseq libraries
    //

    RIBOSEQC_RNA (
        STAR_ALIGN_RNA.out.bam_sorted_aligned,
        ch_bsgenome,
        ch_gtf_Rannot
    )
    ch_versions = ch_versions.mix(RIBOSEQC_RNA.out.versions.first())

    //
    // align Ribo reads with STAR
    //

    STAR_ALIGN_RIBO (
        ch_reads_by_type.ribo,
        ch_genome_index.map { file -> [ [:], file ] },
        ch_gtf.map { file -> [ [:], file ] },
        "$projectDir/assets/NO_FILE", // empty arguments for additional_junctions
        false, // star_ignore_sjdbgtf
        "", // seq_platform
        "" // seq_center
    )
    ch_versions = ch_versions.mix(STAR_ALIGN_RIBO.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN_RIBO.out.log_final.collect{ _meta, log -> log })

    //
    // Quality control with Ribo-seQC for RiboSeq libraries
    //

    RIBOSEQC_RIBO (
        STAR_ALIGN_RIBO.out.bam_sorted_aligned,
        ch_bsgenome,
        ch_gtf_Rannot
    )
    ch_versions = ch_versions.mix(RIBOSEQC_RIBO.out.versions.first())

    //
    // Convert bedGraph to bigWig
    //

    UCSC_BEDGRAPHTOBIGWIG (
        RIBOSEQC_RNA.out.bedgraph
            .mix( RIBOSEQC_RIBO.out.bedgraph )
            .map { _meta, files -> files }
            .flatten()
            .map { file -> [ [ id: file.name.replaceFirst(/\.bedgraph$/, '') ], file ] },
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //

    STAR_ALIGN_RNA.out.bam_sorted_aligned
        .mix( STAR_ALIGN_RIBO.out.bam_sorted_aligned )
        .set { ch_final_bam_aligned }

    SAMTOOLS_INDEX ( ch_final_bam_aligned )
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX.out.versions.first() )

    ch_final_bam_aligned
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }

    BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta.map { file -> [ [:], file ] }.combine(ch_fai).collect() )

    ch_multiqc_files  = ch_multiqc_files.mix( BAM_STATS_SAMTOOLS.out.stats.collect{ _meta, file -> file } )
        .mix( BAM_STATS_SAMTOOLS.out.flagstat.collect{ _meta, file -> file } )
        .mix( BAM_STATS_SAMTOOLS.out.idxstats.collect{ _meta, file -> file } )

    //
    // Merge bams and counts regions files from Ribo-seQC for RNA and Ribo libraries
    //
    counts_regions = RIBOSEQC_RNA.out.counts_regions
        .mix( RIBOSEQC_RIBO.out.counts_regions )

    emit:
    bam              = ch_final_bam_aligned                // channel: [ val(meta), path(bam) ]
    counts_regions   = counts_regions                      // channel: [ val(meta), path(counts_regions) ]
    bam_for_orfquant = RIBOSEQC_RIBO.out.bam_for_orfquant  // channel: [ val(meta), path(bam_for_orfquant) ]
    multiqc_files    = ch_multiqc_files                    // channel: [ logs ]
    versions         = ch_versions                         // channel: [ versions.yml ]
}
