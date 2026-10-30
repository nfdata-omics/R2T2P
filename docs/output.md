# nfdata-omics/r2t2p: Output

## Introduction

This document describes all the outputs produced by the pipeline. 

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Genome index generation](#genome-index-generation) - Generation of the genome index with the aligner *STAR*
- [First alignment](#first-alignment) - First RNA-seq alignment with *STAR* in two-pass mode
- [Second alignment](#second-alignment) - Second RNA-seq alignment with *STAR*
- [Transcriptome assembly](#transcriptome-assembly) - Reference-guided transcriptome reconstruction with *StringTie3*, followed by merging between transcriptome annotations
- [Final alignment](#final-alignment) - Final alignment by providing *STAR* with RNA-seq and/or Ribo-seq data
- [Differential transcriptomic analyses](#differential-transcriptomic-analyses) - Differential RNA-seq and/or Ribo-seq analyses with the R packages *DESeq2* and *DEXSeq*
- [ORF detection](#orf-detection) - Isoform-level *de novo* ORF detection with the R package *ORFquant*, followed by ORF-level differential analyses in R 
- [Proteomic searches and differential analyses](#proteomic-searches-and-differential-analyses) - Proteomic database searches and quantification with the computational platform *FragPipe*, followed by proteomic differential analyses at multiple levels
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Genome index generation

<details markdown="1">
<summary>Output files</summary>

- `genome/star/`
  - `chrLength.txt`: Chromosome lengths
  - `chrName.txt`: Chromosome names
  - `chrNameLength.txt`: Combined chromosome info
  - `chrStart.txt`: Chromosome start positions
  - `Genome`: Packed genome sequence
  - `genomeParameters.txt`: Index generation parameters
  - `SA`: Suffix array binary file
  - `SAindex`: Suffix array index

</details>

[Genome index generation](#genome-index-generation) is performed with *STAR 2.7.9a* (https://github.com/alexdobin/STAR/releases?page=2#release-2.7.9a), which is also used for all alignment steps of the pipeline.

### First alignment

<details markdown="1">
<summary>Output files and directories</summary>

- `first_alignment/`
  - One subdirectory per sample, named after the sample and containing the following files:
    - `[sample_name].first_alignment.Aligned.sortedByCoord.out.bam`: Sorted by coordinate BAM file
    - `[sample_name].first_alignment.Aligned.sortedByCoord.out.bam.bai`: Index file for the BAM file
    - `[sample_name].first_alignment.flagstat`, produced by *samtools flagstat* (https://www.htslib.org/doc/samtools-flagstat.html): Counts for each of 13 categories based primarily on bit flags in the FLAG field. Information on the meaning of the flags is given in the SAM specification document (https://samtools.github.io/hts-specs/SAMv1.pdf).
    - `[sample_name].first_alignment.idxstats`, produced by *samtools idxstats* (https://www.htslib.org/doc/samtools-idxstats.html): Alignment summary statistics, with each line consisting of reference sequence name, sequence length, number of mapped read-segments, and number of unmapped read-segments.
    - `[sample_name].first_alignment.Log.final.out`: Summary mapping statistics after mapping job is complete (very useful for quality control). 
    - `[sample_name].first_alignment.Log.out`: Main log file with a lot of detailed information about the run (useful for troubleshooting and debugging).
    - `[sample_name].first_alignment.Log.progress.out`: This file reports job progress statistics.
    - `[sample_name].first_alignment.SJ.out.tab`: This file contains high confidence collapsed splice junctions in tab-delimited format.
    - `[sample_name].first_alignment.stats`, produced by *samtools stats* (https://www.htslib.org/doc/samtools-stats.html): Comprehensive statistics from the alignment file, including read lengths and mapping qualities.

</details>

The [first alignment](#first-alignment) is performed by providing the aligner *STAR* with RNA-seq data. The aligner is used in two-pass mode. 

### Second alignment

<details markdown="1">
<summary>Output files and directories</summary>

`second_alignment/`
  - One subdirectory per sample, named after the sample and containing the following files:
    - The same types of files as those contained in the sample-specific subdirectories of [First alignment](#first-alignment), with "first_alignment" replaced by "second_alignment" in the file names. 
    - `[sample_name].second_alignment.Aligned.sortedByCoord.out.bam_counts_regions`
    - bigWig files with coverage values on the minus and plus strand, using all reads or only uniquely mapped reads (`[sample_name].second_alignment.Aligned.sortedByCoord.out.bam_coverage_minus.bigWig`,
      `[sample_name].second_alignment.Aligned.sortedByCoord.out.bam_coverage_plus.bigWig`, `[sample_name].second_alignment.Aligned.sortedByCoord.out.bam_coverage_uniq_minus.bigWig`, `[sample_name].second_alignment.Aligned.sortedByCoord.out.bam_coverage_uniq_plus.bigWig`)
    - `[sample_name].second_alignment.Aligned.sortedByCoord.out.bam_readsRle.RData`
    - `[sample_name].second_alignment.Aligned.sortedByCoord.out.bam_strandedness`

</details>

The [second alignment](#second-alignment) is performed by providing the aligner *STAR* with RNA-seq data. Splice junction discovery is disabled.

### Transcriptome assembly

<details markdown="1">
<summary>Output files and directories</summary>

`transcriptome_assembly/`
 - `combined_annotations_Rannot`: Rannot file of the output transcriptome annotation, which is obtained by combining annotations
 - `combined_annotations.gtf`: GTF file of the output transcriptome annotation, which is obtained by combining annotations
 - *GFFCompare* output files of the comparison with the reference annotation (`comp_denovo.annotated.gtf`, `comp_denovo.loci`, `comp_denovo.merged_bams.annotated_ok.gtf.refmap`, 
   `comp_denovo.merged_bams.annotated_ok.gtf.tmap`, `comp_denovo.stats`, and `comp_denovo.tracking`). The .refmap and .tmap files are present only when the user-provided annotation 
   is not specified and used.
 - `comp_denovo.annotated_ok.gtf`: Filtered version of the GTF file `comp_denovo.annotated.gtf` (retaining only regions with known strand)
 - `merged_bams/`
   - `merged_bams.coverage.gtf`: *StringTie3* output GTF file with fully covered transcripts that match the reference annotation
   - `merged_bams.gene.abundance.txt`: *StringTie3* output text file with gene abundances
   - `merged_bams.transcripts.gtf`: *StringTie3* output GTF file containing all the assembled transcripts
   - `merged_bams.annotated_ok.gtf`: Filtered version of the GTF file `merged_bams.transcripts.gtf` (retaining only regions with known strand)
   - `merged_bams.ballgown/`: *StringTie3* output table files (.ctab files), which are required as input to *Ballgown*
 
 When the user-provided annotation is specified and used, the following files and directories are also present in `transcriptome_assembly`:
 - *GFFCompare* output files of the comparison with the reference annotation (`comp_denovo.comp_strg_user.merged_with_strg.gtf.refmap` and `comp_denovo.comp_strg_user.merged_with_strg.gtf.tmap`)
 - `merge_with_user_provided/`
   - *GFFCompare* output files of the comparison with the *StringTie3* annotation (`comp_strg_user.annotated.gtf`, `comp_strg_user.[user_gtf_filename].refmap`, `comp_strg_user.[user_gtf_filename].tmap`, 
   `comp_strg_user.loci`, `comp_strg_user.stats`, `comp_strg_user.tracking`)
   - `comp_strg_user.annotated_ok.gtf`: Filtered version of the GTF file `comp_strg_user.annotated.gtf` (retaining only regions with known strand)
   - `comp_strg_user.merged_with_strg.gtf`: GTF file which is obtained by combining the *StringTie3* and the user-provided transcriptome annotations
   - `gene_ID_pairs.tsv`: Text file with details on which original *StringTie3* and user gene IDs were combined when defining the unified gene identifiers.

</details>
  
[Transcriptome assembly](#transcriptome-assembly) is performed with *StringTie v3.0.0*, a major update of the widely used transcriptome assembler *StringTie* (https://github.com/gpertea/stringtie). 
Following transcriptome reconstruction, annotations are combined by using *GFFCompare* (https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) as well as custom R scripts. 

### Final alignment

<details markdown="1">
<summary>Output files and directories</summary>

`final_alignment/`
 - One subdirectory per sample, named after the sample and containing the same types of files as those contained in the sample-specific subdirectories of [Second alignment](#second-alignment), with some filename 
   changes: “.second_alignment” is replaced by “_RNA.final_alignment” or “_Ribo.final_alignment” for the output files of *samtools flagstat*, *samtools idxstats* and *samtools stats*, and by “_RNA” or “_Ribo” for 
   all other files. 
   Only for Ribo-seq alignments the following files are also present in each sample-specific subdirectory:
   - `[sample_name]_Ribo.Aligned.sortedByCoord.out.bam_for_ORFquant`: *RiboseQC* output file which is provided as input to the R package *ORFquant*
   - `[sample_name]_Ribo.Aligned.sortedByCoord.out.bam_P_sites_calcs`
   - bigWig files for P-sites on the minus and plus strand, using all reads or only uniquely mapped reads (`[sample_name]_Ribo.Aligned.sortedByCoord.out.bam_P_sites_minus.bigWig`, `[sample_name]_Ribo.Aligned.sortedByCoord.out.bam_P_sites_plus.bigWig`,
     `[sample_name]_Ribo.Aligned.sortedByCoord.out.bam_P_sites_uniq_minus.bigWig`, `[sample_name]_Ribo.Aligned.sortedByCoord.out.bam_P_sites_uniq_plus.bigWig`)
   - `[sample_name]_Ribo.Aligned.sortedByCoord.out.bam_results_RiboseQC_all`: R file including all the *RiboseQC* results
   - `[sample_name]_Ribo.Aligned.sortedByCoord.out.bam_results_RiboseQC`: R file including some *RiboseQC* results 
   - `[sample_name]_Ribo.ReadsPerGene.out.tab`: *STAR* output tab-separated file with read counts per gene
   
</details>

The [final alignment](#final-alignment) is performed by providing the aligner *STAR* with RNA-seq and/or Ribo-seq data.

### Differential transcriptomic analyses

<details markdown="1">
<summary>Output files and directories</summary>

`differential_expression/`
 - One subdirectory per comparison, named after the comparison (i.e., [comparison] = [condition]\_vs\_[baseline]) and containing the following files:
   - `[comparison]_countmats.RData`: Count matrices
   - `[comparison]_multiDE_annotfiles.RData`
   - `[comparison]_multiDE_input.RData`: RData file which is provided as input to multiDE
   - Results of multiDE, RF regression and GO enrichment analysis (`[comparison]_multiDE_results.RData`, `[comparison]_multiDE_results_RF.RData`, `[comparison]_multiDE_results_summary_allgenes_all.tsv`,
    `[comparison]_multiDE_results_summary_allgenes.RData`, `[comparison]_multiDE_results_summary_allgenes_uniq.tsv`, `[comparison]_multiDE_results_summary_all.tsv`, `[comparison]_multiDE_results_summary_GO.RData`,
    `[comparison]_multiDE_results_summary.RData`, `[comparison]_multiDE_results_summary_subgene_all.tsv`, `[comparison]_multiDE_results_summary_subgene.RData`, `[comparison]_multiDE_results_summary_subgene_uniq.tsv`,    `[comparison]_multiDE_results_summary_uniq.tsv`)
  
</details>

[Differential transcriptomic analyses](#differential-transcriptomic-analyses) are performed at gene level and at subgene level (exon bins, junctions, introns) by using the R packages *DESeq2* (https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and *DEXSeq* (https://bioconductor.org/packages//release/bioc/html/DEXSeq.html).

### ORF detection

<details markdown="1">
<summary>Output files and directories</summary>

`ORFquant/`
 - `combined_forORFquant.RData_Detected_ORFs.gtf`: GTF file containing coordinates of the detected ORFs
 - `combined_forORFquant.RData_final_ORFquant_results`: RData object file containing the final ORFquant results
 - `combined_forORFquant.RData_Protein_sequences.fasta`: FASTA file containing the set of translated proteins
 - `combined_forORFquant.RData_tmp_ORFquant_results`: RData object file containing the entire set of results for each genomic region
 - A subdirectory for each comparison, named after the comparison (i.e., [comparison] = [condition]\_vs\_[baseline]) and containing the 
   following files:
   - `[comparison]_ORFDEX.RData`: RData file with results of ORF-level differential analyses
   - `[comparison]_ORFquant_counts.RData`: RData file with counts
   - `[comparison]_ORFquant_feats.RData`: RData file with features
   
</details>

[ORF detection](#orf-detection) is performed by providing the R package *ORFquant* (https://github.com/lcalviell/ORFquant) with P-sites positions, which are extracted from Ribo-seq data using *RiboseQC* (https://github.com/ohlerlab/RiboseQC).

### Proteomic searches and differential analyses

<details markdown="1">
<summary>Output files and directories</summary>

`fragpipe/`
 - The following files and directories are present for each proteomic search, with [database] replaced by the specific database that is used for the search (i.e., [database] = 'ref_database', 'custom_database', or 'combined_database'):
   - `[database].fasta`: FASTA file of the protein database (only target protein sequences, without contaminants and decoys)
   - A subdirectory `[database]_search`: Folder containing proteomic search results
   - `[database]_search_peptide_annotation.RData`: RData file with peptide-to-proteins mapping information
   - `[database]_search_peptide_coords.gtf`: GTF file of peptides
   - `[database]_search_peptide_intensities.tsv`: TSV file with peptide intensities
   - bigWig files with peptide log2 fold changes for each comparison (`[database]_search_peptide_log2FC_[comparison]_minus.bw` and `[database]_search_peptide_log2FC_[comparison]_plus.bw`)
   - `[database]_search_proteomics_multiDE.RData`: RData file with results of proteomic differential analyses with the R package *limpa*
   - Additional RData files with results of proteomic analyses (`[database]_search_peptide_results_annotated.RData` and `[database]_search_proteomics_all_intensities.RData`)
   - `[database].workflow`: FragPipe workflow file 
   - `decoys-contam-[database].fasta.fas`: FASTA file of the protein database after adding decoys and contaminants (target protein sequences, contaminants, and decoys)
   - *MSFragger* peptide index files (`decoys-contam-[database].fasta.fas.1.pepindex` and `decoys-contam-[database].fasta.fas.2.pepindex`)

</details>

[Proteomic searches and differential analyses](#proteomic-searches-and-differential-analyses) are performed by providing experimental data and protein databases to the computational platform *FragPipe*, which combines multiple well-established tools within a single framework. Various combinations of data analysis tools and different search settings can be used, depending on the data type (DIA, DDA-LFQ, or DDA-TMT) as well as on user choices. 

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report. Results generated by *MultiQC* collate pipeline QC from supported tools. Software versions are reported in the *MultiQC* output for future traceability. For more information about how to use *MultiQC* reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_[YYYY-MM-DD_HH-MM-SS].html`: HTML report including metrics about a pipeline run. It is organized into three sections: Summary, Resources, and Tasks. 
  - `execution_timeline_[YYYY-MM-DD_HH-MM-SS].html`: HTML report showing a timeline of the tasks executed in a pipeline run
  - `execution_trace_[YYYY-MM-DD_HH-MM-SS].txt`: Text file containing information about each task executed in a pipeline run, including name, status, submission time, and duration
  - `params_[YYYY-MM-DD-HH-MM-SS].json`: JSON file containing parameters used by the pipeline run
  - `pipeline_dag_[YYYY-MM-DD_HH-MM-SS].html`: HTML file containing a representation of the pipeline as a directed acyclic graph (DAG)
  - `pipeline_report.html`: HTML report summarizing a pipeline execution, including the run name, the command used to launch it, and pipeline configuration details
  - `pipeline_report.txt`: Text report summarizing a pipeline execution, including the run name, the command used to launch it, and pipeline configuration details (the same information as in pipeline_report.html, but in a text file)
  - `r2t2p_software_mqc_versions.yml`: YAML file listing software versions used in different pipeline steps

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
