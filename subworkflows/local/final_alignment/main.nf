//
// Final alignment against the new transcriptome
//

include { STAR_ALIGN as STAR_ALIGN_RNA                    } from '../../../modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_RIBO                   } from '../../../modules/nf-core/star/align/main'
include { BAM_SORT_STATS_SAMTOOLS as RNA_BAM_SORT_STATS   } from '../../../subworkflows/nf-core/bam_sort_stats_samtools'
include { BAM_SORT_STATS_SAMTOOLS as RIBO_BAM_SORT_STATS  } from '../../../subworkflows/nf-core/bam_sort_stats_samtools'

workflow FINAL_ALIGNMENT {
    take:
    ch_reads        // channel: [ val(meta), [ fastq ] ]
    ch_genome_index // channel: genome_index
    ch_gtf          // channel: gtf
    ch_fasta        // channel: fasta

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_reads
        .branch { meta, _fastq ->
            rna:  meta.library_type == "RNA"
            ribo: meta.library_type == "Ribo"
        }
    .set { ch_reads_by_type }

    // align RNA reads with STAR
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

// $star_path --genomeDir $genome_index --readFilesIn $full_fasta --runThreadN $SLURM_CPUS_PER_TASK --outFilterMismatchNmax 4
//     --outFilterMultimapNmax 250 --chimScoreSeparation 10 --chimScoreMin 20 --chimSegmentMin 15
//     --outSAMattributes NH HI AS nM NM MD XS --alignSJoverhangMin 500 --outFileNamePrefix $name_exp"_"
//      --outReadsUnmapped Fastx --outWigNorm RPM --outSAMstrandField intronMotif --outSAMmultNmax 1
//      --outMultimapperOrder Old_2.4 --limitBAMsortRAM 19000000000 --alignEndsType "Local"
//      --readFilesCommand zcat --sjdbGTFfile $stringtie_gtf_file --outSAMtype BAM SortedByCoordinate --limitSjdbInsertNsj 5000000

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    RNA_BAM_SORT_STATS ( STAR_ALIGN_RNA.out.bam, ch_fasta.map { file -> [ [:], file ] } )
    ch_versions = ch_versions.mix(RNA_BAM_SORT_STATS.out.versions)
    ch_multiqc_files  = ch_multiqc_files.mix( RNA_BAM_SORT_STATS.out.stats.collect{ _meta, log -> log } )
        .mix( RNA_BAM_SORT_STATS.out.flagstat.collect{ _meta, log -> log } )
        .mix( RNA_BAM_SORT_STATS.out.idxstats.collect{ _meta, log -> log } )

// bamfull="`readlink -f $name_exp"_Aligned.sortedByCoord.out.bam"`"
// $samtools_path index $bamfull
// Rscript $6"/f_run_newriboseQC_cnt.R" $bamfull ../../transcriptome_assembly/comp_denovo.annotated_ok.gtf_stringtie_Rannot $6
// find . -name "*bedgraph" | sort | xargs -n 1 sh -c 'sbatch -t 01:00:00 --mem=40000  ../../f_bedgraph_tobw.bash $0'


    // (Final alignment steps would go here)

    emit:
    gff_stats = ch_multiqc_files                    // channel: [ stats ]
    versions  = ch_versions                         // channel: [ versions.yml ]
}
