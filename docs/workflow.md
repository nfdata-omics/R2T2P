# nfdata-omics/r2t2p: Workflow Rationale

## Introduction

R2T2P is designed to integrate RNA-seq, Ribo-seq, isoform-level ORF discovery, and proteomics into a
single reproducible workflow. The central rationale is that no single data type fully resolves the expressed
and translated transcriptome. RNA-seq provides evidence for transcript structures and gene-level abundance,
Ribo-seq provides evidence for translation, and LC-MS/MS proteomics can validate translated products at
peptide level.

The workflow is organized into three main biological analysis modules:

- The RNA module reconstructs and refines transcript annotations from short-read RNA-seq data.
- The Translation module identifies translated regions on annotated and newly reconstructed transcripts using
  Ribo-seq evidence.
- The Protein module searches mass-spectrometry data against custom protein databases derived from reference
  CDSs and ORFquant-predicted ORFs.

R2T2P also performs quality-control checks and differential analyses across the supported data layers. The
outputs are collected into reports that summarize the full workflow, including read QC, mapping statistics,
annotation-comparison metrics, Ribo-seq processing statistics, proteomics search information, and software
versions.

## Pipeline overview

R2T2P builds progressively from reference-guided alignment to annotation refinement, then uses the refined
annotation to quantify transcription, ribosome occupancy, translated ORFs, and protein evidence. This staged
design keeps the curated reference annotation as a stable backbone while allowing sample-specific or
user-provided transcript models to be added when they are supported by the data.

The major workflow stages are:

- Prepare reference files, input reads, and quality-control reports.
- Align RNA-seq reads to the genome to identify splice junctions and support transcript discovery.
- Assemble and merge transcript annotations using RNA-seq evidence, the reference annotation, and optionally a
  user-provided GTF.
- Re-align RNA-seq and Ribo-seq reads against the refined annotation.
- Quantify transcriptomic and ribosome-profiling signal for differential analyses.
- Detect translated regions with ORFquant from Ribo-seq evidence.
- Build protein databases from annotated CDSs and predicted ORFs for optional proteomic searches.
- Summarize quality-control metrics, software versions, and workflow provenance with MultiQC.

## RNA module

The RNA module reconstructs a transcriptome annotation from short-read RNA-seq data. The goal is to improve
transcript-level interpretation without discarding the reliability of a curated reference annotation.

Before transcriptome assembly, two STAR genome alignments are performed. In the first alignment, RNA-seq reads
are aligned with the supplied reference annotation to guide splice-junction detection while still allowing
STAR to discover exon-exon junctions that are absent from the annotation. Novel junctions detected across RNA-seq
samples are then collected and filtered.

The second RNA-seq alignment uses the reference annotation together with the junctions identified in the first
alignment. This reduces the chance that downstream assembly is driven by unsupported splice junctions while
still preserving sample-specific junction evidence. The resulting alignments are sorted, indexed, and assessed
with standard mapping statistics.

The RNA module then uses StringTie to assemble transcripts from the RNA-seq alignments, with the reference
annotation supplied to guide transcript reconstruction. This enables detection of novel isoforms in complex
transcriptomes while keeping known transcript models available during assembly.

In the final stage of the RNA module, the assembled annotation is compared to the reference annotation with
GFFCompare. GFFCompare class codes are used to describe assembled transcripts relative to reference transcripts
and to report performance metrics at exon, intron, intron-chain, and transcript levels. Custom R scripts then
merge the assembled and reference annotations. Transcript biotypes and GFFCompare classes are harmonized during
this step.

StringTie-assembled transcripts that are already represented in the reference annotation are excluded from the
novel transcript set to avoid redundancy. This includes exact intron-chain matches and transcripts that are
fully contained within reference transcripts and intron-compatible with them. If the user provides an additional
GTF annotation, the pipeline can first compare and merge that annotation with the assembled transcriptome, then
merge the result with the reference annotation.

The RNA module produces an augmented transcriptome annotation containing reference transcripts and supported
novel transcripts. It also creates analysis-ready annotation objects and feature sets, including exons, CDSs,
UTRs, introns, junctions, intergenic regions, exonic bins, gene-level CDS coordinates, transcript-level CDS
coordinates, start codons, and stop codons.

RNA module quality checks include:

- STAR mapping statistics for both RNA-seq alignment steps.
- samtools alignment summaries for sorted and indexed BAM files.
- GFFCompare metrics for transcriptome assembly and annotation comparison.

## Translation module

The Translation module identifies translated regions at isoform level. It uses RNA-seq and Ribo-seq reads
together with the augmented transcriptome annotation generated by the RNA module.

First, RNA-seq and Ribo-seq reads are aligned to the genome with STAR using the exon-exon junctions and
transcript models from the augmented annotation. Running both assays against the same refined annotation places
transcriptional and translational evidence in a consistent coordinate system.

RiboseQC is then used to process the final alignment files. It generates read counts and statistics across
annotation features such as genes, coding regions, untranslated regions, and other derived regions. These
outputs provide the bridge between genome alignments and downstream R analyses, including differential analyses
and ORFquant input generation.

ORFquant uses Ribo-seq alignments and RiboseQC outputs to identify actively translated regions on transcripts.
The isoform-aware design is important because alternative isoforms of the same gene can have different coding
potential, untranslated regions, or translation levels. This makes ORFquant suitable for studying unannotated
translated regions in complex transcriptomes.

R2T2P also supports differential analyses from RNA-seq and Ribo-seq quantifications. Gene-level analyses are
performed from read-count summaries, and ORF-level differential analyses are performed from ORFquant results
when contrasts are provided. Interpreting RNA-seq and Ribo-seq together helps distinguish changes in RNA
abundance from changes in ribosome occupancy or translated-region usage.

The Translation module also produces protein FASTA databases for proteomics:

- A database of annotated protein sequences extracted from reference CDS coordinates.
- A database of protein sequences corresponding to ORFquant-detected ORFs.
- A combined database containing both annotated proteins and ORFquant-derived protein sequences.

The annotated database supports searches against known proteins, the ORFquant database focuses searches on
candidate novel translation products, and the combined database enables a broader search that includes both
evidence sources.

Translation module quality checks include:

- STAR mapping statistics for final RNA-seq and Ribo-seq alignments.
- Read counts, read lengths, and annotation-feature statistics from RiboseQC.
- Numbers, lengths, annotation status, and quantification estimates of ORFquant-detected ORFs.

## Protein module

The Protein module extends the workflow from inferred translation to peptide-level evidence. It can run
FragPipe searches when mass-spectrometry inputs, a FragPipe manifest, a workflow file, and the required external
tools are provided.

Before searching, decoys and contaminants are added to the protein databases with Philosopher. This allows
peptide-spectrum matches, peptides, and proteins to be assessed with standard false-discovery-rate procedures.

FragPipe is then run with the user-provided workflow file. This design allows R2T2P to support different
proteomics strategies through FragPipe configuration, including DDA, TMT, label-free quantification, and DIA
workflows when the corresponding FragPipe settings and external tools are supplied. Depending on the workflow
file, FragPipe can combine tools such as MSFragger, MSFragger-DIA, MSBooster, Percolator, IonQuant,
TMT-Integrator, EasyPQP, DIA-Umpire, DIA-NN, and Philosopher.

Search results can be interpreted at peptide, protein-group, and gene levels. Comparing searches against the
annotated, ORFquant-derived, and combined databases helps evaluate whether candidate novel ORFs are supported by
proteomic evidence and whether spectra are assigned differently when novel translation products are available
in the search space.

Protein module quality checks include:

- Numbers of detected peptide-spectrum matches across samples.
- Numbers of detected and quantified peptides, protein groups, and genes.
- FDR filtering settings and thresholds for peptide-spectrum matches, peptides, ions, and proteins.
- FragPipe log files documenting execution of the proteomics analysis steps.
