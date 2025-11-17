//
// Final alignment against the new transcriptome
//

include { STAR_ALIGN as STAR_ALIGN_RNA                    } from '../../../modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_RIBO                   } from '../../../modules/nf-core/star/align/main'
include { BAM_SORT_STATS_SAMTOOLS as RNA_BAM_SORT_STATS   } from '../../../subworkflows/nf-core/bam_sort_stats_samtools'
include { BAM_SORT_STATS_SAMTOOLS as RIBO_BAM_SORT_STATS  } from '../../../subworkflows/nf-core/bam_sort_stats_samtools'
include { RIBOSEQC as RIBOSEQC_RNA                        } from '../../../modules/local/riboseqc/main'
include { RIBOSEQC as RIBOSEQC_RIBO                       } from '../../../modules/local/riboseqc/main'
include { UCSC_BEDGRAPHTOBIGWIG                           } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow FINAL_ALIGNMENT {
    take:
    ch_reads        // channel: [ val(meta), [ fastq ] ]
    ch_genome_index // channel: genome_index
    ch_gtf          // channel: gtf
    ch_fasta        // channel: fasta
    ch_chrom_sizes  // channel: chrom_sizes
    ch_bsgenome     // channel: bsgenome
    ch_gtf_Rannot   // channel: gtf_Rannot

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_reads
        .branch { meta, _fastq ->
            rna:  meta.library_type == "RNA"
            ribo: meta.library_type == "Ribo"
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
    // Convert bedGraph to bigWig
    //
    UCSC_BEDGRAPHTOBIGWIG (
        RIBOSEQC_RNA.out.bedgraph
            .map { _meta, files -> files}
            .flatten()
            .map { file -> [ [ id: file.name.replaceFirst(/\.bedgraph$/, '') ], file ] },
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())

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
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    STAR_ALIGN_RNA.out.bam_sorted_aligned
        .mix( STAR_ALIGN_RIBO.out.bam_sorted_aligned )
        .set { ch_final_bam_aligned }
    RNA_BAM_SORT_STATS ( ch_final_bam_aligned, ch_fasta.map { file -> [ [:], file ] } )
    ch_versions = ch_versions.mix(RNA_BAM_SORT_STATS.out.versions)
    ch_multiqc_files  = ch_multiqc_files.mix( RNA_BAM_SORT_STATS.out.stats.collect{ _meta, log -> log } )
        .mix( RNA_BAM_SORT_STATS.out.flagstat.collect{ _meta, log -> log } )
        .mix( RNA_BAM_SORT_STATS.out.idxstats.collect{ _meta, log -> log } )

    //
    // Merge counts regions files from Ribo-seQC for RNA and Ribo libraries
    //
    counts_regions = RIBOSEQC_RNA.out.counts_regions
        .mix( RIBOSEQC_RIBO.out.counts_regions )

    emit:
    counts_regions = counts_regions                     // channel: [ val(meta), path(counts_regions) ]
    multiqc_files = ch_multiqc_files                    // channel: [ logs ]
    versions      = ch_versions                         // channel: [ versions.yml ]
}
