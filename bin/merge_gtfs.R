#!/usr/bin/Rscript

# This script merges reference and StringTie GTF annotations, processes transcript and gene features,
# and generates a combined annotation object for downstream analysis. It also exports a merged GTF file
# and saves an R annotation object for use in RiboseQC or similar workflows.

suppressPackageStartupMessages(library("RiboseQC"))

# Parse arguments: reference annotation RData, StringTie GTF, and a logical flag for merging
args <- commandArgs(trailingOnly = TRUE)

# Load reference annotation and create TxDb from StringTie GTF
GTF_annotation <- get(load(args[1]))
annotation_ok <- GTF_annotation
annotation_str <- makeTxDbFromGFF(file = args[2])

# Extract exons and genes from both reference and StringTie annotations
exs_annot <- GTF_annotation$exons_txs
exs_str <- exonsBy(annotation_str, by = "tx", use.names = TRUE)
genes_annot <- GTF_annotation$genes
genes_str <- GenomicFeatures::genes(annotation_str)

# Map transcripts to genes for both annotations
txs_gene_map <- GTF_annotation$trann
transcripts_db_str <- transcripts(annotation_str)
txs_gene_map_str <- mapIds(
  transcripts_db_str$tx_name,
  x = annotation_str,
  column = "GENEID",
  keytype = "TXNAME"
)
txs_gene_map_str <- data.frame(
  gene_id = txs_gene_map_str,
  transcript_id = names(txs_gene_map_str),
  stringsAsFactors = FALSE
)
rownames(txs_gene_map_str) <- NULL

# Process StringTie annotation: extract transcript/gene info and assign biotypes/class codes
trann <- unique(mcols(import.gff2(
  args[2],
  colnames = c(
    "gene_id", "gene_biotype", "gene_type", "gene_name", "gene_symbol",
    "transcript_id", "transcript_biotype", "transcript_type",
    "cmp_ref", "class_code"
  )
)))
trann <- trann[!is.na(trann$class_code), ]
trann <- data.frame(unique(trann), stringsAsFactors = FALSE)

# Handle missing biotype/name/type fields and assign default values if necessary
# (Several if statements below)

# Assign class codes and biotypes to StringTie transcripts
colnames(trann) <- c(
  "gene_id", "gene_biotype", "gene_name",
  "transcript_id", "transcript_biotype", "nearest_ref", "class_code"
)
str_ann <- DataFrame(trann)
trann <- GTF_annotation$trann

# Adjust class codes for intronic transcripts with different strand than host
# Split StringTie annotation into duplicate, new, and isoform categories

# Map class codes to descriptive biotype labels
code_annot <- cbind(
  c("=", "c", "k", "r", "u", "x", "s", "e", "i", "j", "o", "p", "y", "m", "n", "xi"),
  paste("novel_", c(
    "same", "contained", "extended", "repeat", "intergenic", "antisense", "antisense_spl",
    "pre_mRNA", "intronic", "alt_isof", "exon_overl", "run_on", "enclosing",
    "retain_intr", "retain_partintr", "antisense_intr"
  ), sep = "")
)

# Merge StringTie annotation with reference, depending on the logical flag
if (as.logical(args[3]) == FALSE) {
  trann <- rbind(trann, str_ann[!str_ann$class_code %in% c("c", "="), colnames(trann)])
} else if (as.logical(args[3]) == TRUE) {
  trann <- rbind(trann, str_ann[, colnames(trann)])
}

# Fill missing gene names with gene IDs
trann$gene_name[is.na(trann$gene_name)] <- trann$gene_id[is.na(trann$gene_name)]

all_names <- names(c(exs_str, exs_annot))
missing <- setdiff(trann$transcript_id, all_names)
print(all_names[1:10])
print(trann$transcript_id[1:10])

# Combine exons and CDS from both annotations, assign metadata columns
exs_gtf <- unlist(c(exs_str, exs_annot)[trann$transcript_id])
mcols(exs_gtf) <- NULL
exs_gtf$transcript_id <- names(exs_gtf)
exs_gtf$gene_id <- trann[match(exs_gtf$transcript_id, trann$transcript_id), "gene_id"]
exs_gtf$gene_biotype <- trann[match(exs_gtf$gene_id, trann$gene_id), "gene_biotype"]
exs_gtf$gene_name <- trann[match(exs_gtf$gene_id, trann$gene_id), "gene_name"]
exs_gtf$transcript_biotype <- trann[match(exs_gtf$transcript_id, trann$transcript_id), "transcript_biotype"]
exs_gtf$type <- "exon"

cds_gtf <- unlist(GTF_annotation$cds_txs)
mcols(cds_gtf) <- NULL
cds_gtf$transcript_id <- names(cds_gtf)
cds_gtf$gene_id <- trann[match(cds_gtf$transcript_id, trann$transcript_id), "gene_id"]
cds_gtf$gene_biotype <- trann[match(cds_gtf$gene_id, trann$gene_id), "gene_biotype"]
cds_gtf$gene_name <- trann[match(cds_gtf$gene_id, trann$gene_id), "gene_name"]
cds_gtf$transcript_biotype <- trann[match(cds_gtf$transcript_id, trann$transcript_id), "transcript_biotype"]
cds_gtf$type <- "CDS"

# Merge and sort all features, check for strand consistency, and filter problematic genes
all <- sort(c(exs_gtf, cds_gtf))
all$source <- paste0(args[2], "_merged_with_ref")
names(all) <- NULL

all2 <- split(all, all$gene_id)
stra <- strand(all2)
genes_ok <- names(which(elementNROWS(runValue(stra)) == 1))
if (sum(elementNROWS(runValue(stra)) > 1) > 0) stop("Some genes have incosistent strand info")
all_ok <- all2[genes_ok]

all <- sort(unlist(all_ok))
names(all) <- NULL

# Export merged annotation as GTF file
gtf_file <- paste0(args[2], "_stringtie.gtf")
export.gff2(object = all, con = gtf_file)

# Create TxDb from merged GTF and extract gene, exon, CDS, UTR, and intron features
seqinfotwob <- GTF_annotation$seqinfo
annotation <- makeTxDbFromGFF(file = gtf_file, format = "gtf", chrominfo = seqinfotwob)

genes <- GenomicFeatures::genes(annotation)
exons_ge <- exonsBy(annotation, by = "gene")
exons_ge <- reduce(exons_ge)

cds_gen <- cdsBy(annotation, "gene")
cds_ge <- reduce(cds_gen)

threeutrs <- reduce(GenomicRanges::setdiff(
  unlist(threeUTRsByTranscript(annotation)),
  unlist(cds_ge),
  ignore.strand = FALSE
))
fiveutrs <- reduce(GenomicRanges::setdiff(
  unlist(fiveUTRsByTranscript(annotation)),
  unlist(cds_ge),
  ignore.strand = FALSE
))
introns <- reduce(GenomicRanges::setdiff(
  unlist(intronsByTranscript(annotation)),
  unlist(exons_ge),
  ignore.strand = FALSE
))
nc_exons <- reduce(GenomicRanges::setdiff(
  unlist(exons_ge),
  reduce(c(unlist(cds_ge), fiveutrs, threeutrs)),
  ignore.strand = FALSE
))

# Assign gene IDs to UTRs, introns, and non-coding exons by overlap
# (findOverlaps and assignment code)

# Define intergenic regions as gaps between genes
intergenicRegions <- genes
strand(intergenicRegions) <- "*"
intergenicRegions <- gaps(reduce(intergenicRegions))
intergenicRegions <- intergenicRegions[strand(intergenicRegions) == "*"]

# Extract transcript, CDS, exon, and intron features by gene and transcript
cds_tx <- cdsBy(annotation, "tx", use.names = TRUE)
txs_gene <- transcriptsBy(annotation, by = "gene")
genes_red <- reduce(sort(GenomicFeatures::genes(annotation)))
exons_tx <- exonsBy(annotation, "tx", use.names = TRUE)
transcripts_db <- transcripts(annotation)
intron_names_tx <- intronsByTranscript(annotation, use.names = TRUE)

# Generate exonic bins, including multi-gene overlaps
nsns <- exonicParts(annotation, linked.to.single.gene.only = FALSE)

# Map CDS coordinates to transcript coordinates for ORF boundary analysis
exsss_cds <- exons_tx[names(cds_tx)]
chunks <- seq(1, length(cds_tx), by = 20000)
if (chunks[length(chunks)] < length(cds_tx)) chunks <- c(chunks, length(cds_tx))
mapp <- GRangesList()
for (i in 1:(length(chunks) - 1)) {
  if (i != (length(chunks) - 1)) {
    mapp <- suppressWarnings(c(
      mapp,
      pmapToTranscripts(
        cds_tx[chunks[i]:(chunks[i + 1] - 1)],
        transcripts = exsss_cds[chunks[i]:(chunks[i + 1] - 1)]
      )
    ))
  }
  if (i == (length(chunks) - 1)) {
    mapp <- suppressWarnings(c(
      mapp,
      pmapToTranscripts(
        cds_tx[chunks[i]:(chunks[i + 1])],
        transcripts = exsss_cds[chunks[i]:(chunks[i + 1])]
      )
    ))
  }
}
cds_txscoords <- unlist(mapp)

# Extract biotypes and IDs from merged GTF for downstream annotation
trann <- unique(mcols(import.gff2(
  gtf_file,
  colnames = c(
    "gene_id", "gene_biotype", "gene_type", "gene_name", "gene_symbol",
    "transcript_id", "transcript_biotype", "transcript_type"
  )
)))
trann <- trann[!is.na(trann$transcript_id), ]
trann <- data.frame(unique(trann), stringsAsFactors = FALSE)

# Handle missing biotype/name/type fields for merged annotation
# (Several if statements below)
colnames(trann) <- c("gene_id", "gene_biotype", "gene_name", "transcript_id", "transcript_biotype")
trann <- DataFrame(trann)

# Assign transcript and gene IDs to introns
unq_intr <- sort(unique(unlist(intron_names_tx)))
names(unq_intr) <- NULL
all_intr <- unlist(intron_names_tx)

ov <- findOverlaps(unq_intr, all_intr, type = "equal")
ov <- split(subjectHits(ov), queryHits(ov))
a_nam <- CharacterList(lapply(ov, function(x) unique(names(all_intr)[x])))

unq_intr$type <- "J"
unq_intr$tx_name <- a_nam

mat_genes <- match(unq_intr$tx_name, trann$transcript_id)
g <- unlist(apply(
  cbind(1:length(mat_genes), Y = elementNROWS(mat_genes)),
  FUN = function(x) rep(x[1], x[2]),
  MARGIN = 1
))
g2 <- split(trann[unlist(mat_genes), "gene_id"], g)
unq_intr$gene_id <- CharacterList(lapply(g2, unique))

# Separate non-coding RNA and non-coding isoform regions
ncrnas <- nc_exons[!nc_exons %over% genes[trann$gene_id[trann$gene_biotype == "protein_coding"]]]
ncisof <- nc_exons[nc_exons %over% genes[trann$gene_id[trann$gene_biotype == "protein_coding"]]]

# Map CDS coordinates to gene IDs and analyze start/stop codons
cds_txscoords$gene_id <- trann$gene_id[match(as.vector(seqnames(cds_txscoords)), trann$transcript_id)]
cds_cc <- cds_txscoords
strand(cds_cc) <- "*"
sta_cc <- resize(cds_cc, width = 1, "start")
sta_cc <- unlist(pmapFromTranscripts(sta_cc, exons_tx[seqnames(sta_cc)], ignore.strand = FALSE))
sta_cc$gene_id <- trann$gene_id[match(names(sta_cc), trann$transcript_id)]
sta_cc <- sta_cc[sta_cc$hit]
strand(sta_cc) <- structure(as.vector(strand(transcripts_db)), names = transcripts_db$tx_name)[names(sta_cc)]
sta_cc$type <- "start_codon"
mcols(sta_cc) <- mcols(sta_cc)[, c("exon_rank", "type", "gene_id")]

sto_cc <- resize(cds_cc, width = 1, "end")
sto_cc <- shift(sto_cc, -2)
stop_inannot <- GTF_annotation$stop_in_gtf
if (is.na(stop_inannot)) {
  sto_cc <- resize(trim(shift(sto_cc, 3)), width = 1, fix = "end")
}
sto_cc <- unlist(pmapFromTranscripts(sto_cc, exons_tx[seqnames(sto_cc)], ignore.strand = FALSE))
sto_cc <- sto_cc[sto_cc$hit]
sto_cc$gene_id <- trann$gene_id[match(names(sto_cc), trann$transcript_id)]
strand(sto_cc) <- structure(as.vector(strand(transcripts_db)), names = transcripts_db$tx_name)[names(sto_cc)]
sto_cc$type <- "stop_codon"
mcols(sto_cc) <- mcols(sto_cc)[, c("exon_rank", "type", "gene_id")]

# Analyze start/stop codon positions: most upstream/downstream, most frequent
start_stop_cc <- sort(c(sta_cc, sto_cc))
start_stop_cc$transcript_id <- names(start_stop_cc)
start_stop_cc$most_up_downstream <- FALSE
start_stop_cc$most_frequent <- FALSE

# Calculate representative transcripts for each gene based on start/stop codon features
# (by and summary code)

cds_txscoords$reprentative_mostcommon <- as.character(seqnames(cds_txscoords)) %in% unlist(repres_freq)
cds_txscoords$reprentative_boundaries <- as.character(seqnames(cds_txscoords)) %in% unlist(repres_upstr)
cds_txscoords$reprentative_5len <- as.character(seqnames(cds_txscoords)) %in% unlist(repres_len5)

# Assign representative status to unique start/stop codon positions
unq_stst <- start_stop_cc
mcols(unq_stst) <- NULL
unq_stst <- sort(unique(unq_stst))
ov <- findOverlaps(unq_stst, start_stop_cc, type = "equal")
ov <- split(subjectHits(ov), queryHits(ov))
unq_stst$type <- CharacterList(lapply(ov, function(x) unique(start_stop_cc$type[x])))
unq_stst$transcript_id <- CharacterList(lapply(ov, function(x) start_stop_cc$transcript_id[x]))
unq_stst$gene_id <- CharacterList(lapply(ov, function(x) unique(start_stop_cc$gene_id[x])))

unq_stst$reprentative_mostcommon <- sum(!is.na(match(unq_stst$transcript_id, unlist(as(repres_freq, "CharacterList"))))) > 0
unq_stst$reprentative_boundaries <- sum(!is.na(match(unq_stst$transcript_id, unlist(as(repres_upstr, "CharacterList"))))) > 0
unq_stst$reprentative_5len <- sum(!is.na(match(unq_stst$transcript_id, unlist(as(repres_len5, "CharacterList"))))) > 0

# Assemble all processed features into a single annotation list for downstream use
GTF_annotation <- list(
  transcripts_db, txs_gene, GTF_annotation$seqinfo, unq_stst, cds_tx, intron_names_tx,
  cds_gen, exons_tx, nsns, unq_intr, genes, threeutrs, fiveutrs, ncisof, ncrnas,
  introns, intergenicRegions, trann, cds_txscoords, GTF_annotation$genetic_codes,
  GTF_annotation$genome_package, GTF_annotation$stop_in_gtf
)
names(GTF_annotation) <- c(
  "txs", "txs_gene", "seqinfo", "start_stop_codons", "cds_txs", "introns_txs",
  "cds_genes", "exons_txs", "exons_bins", "junctions", "genes", "threeutrs",
  "fiveutrs", "ncIsof", "ncRNAs", "introns", "intergenicRegions", "trann",
  "cds_txs_coords", "genetic_codes", "genome_package", "stop_in_gtf"
)

# Warn about transcripts with missing exon boundaries (e.g. trans-splicing events)
txs_all <- unique(GTF_annotation$trann$transcript_id)
txs_exss <- unique(names(GTF_annotation$exons_txs))
txs_notok <- txs_all[!txs_all %in% txs_exss]
if (length(txs_notok) > 0) {
  set.seed(666)
  cat(
    paste(
      "Warning: ", length(txs_notok), " txs with incorrect/unspecified exon boundaries - e.g. trans-splicing events, examples: ",
      paste(txs_notok[sample(1:length(txs_notok), size = min(3, length(txs_notok)), replace = FALSE)], collapse = ", "),
      " - ", date(), "\n", sep = ""
    )
  )
}

# Save the final annotation object for downstream analysis
rannot_file <- paste0(args[2], "_stringtie_Rannot")
save(GTF_annotation, file = rannot_file)
