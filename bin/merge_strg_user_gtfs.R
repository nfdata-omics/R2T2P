#!/usr/bin/Rscript

# Load required package
suppressPackageStartupMessages(library("RiboseQC"))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Load GTF annotations for reference and user-supplied files
GTF_annotation <- makeTxDbFromGFF(file = args[1])
annotation_str <- makeTxDbFromGFF(file = args[2])

# Extract exons and genes from both annotations
exs_annot <- exonsBy(GTF_annotation, by = "tx", use.names = TRUE)
exs_str <- exonsBy(annotation_str, by = "tx", use.names = TRUE)
genes_annot <- GenomicFeatures::genes(GTF_annotation)
genes_str <- GenomicFeatures::genes(annotation_str)

# Build transcript-to-gene mapping for reference annotation
unlist_txs_by_gene <- unlist(GenomicFeatures::transcriptsBy(GTF_annotation, by = "gene"))
gid <- names(unlist_txs_by_gene)
gbt <- rep("novel", length(gid))
gnm <- gid
tid <- unlist_txs_by_gene$tx_name
tbt <- rep("novel", length(tid))

txs_gene_map <- DataFrame(
  gene_id = gid,
  gene_biotype = gbt,
  gene_name = gnm,
  transcript_id = tid,
  transcript_biotype = tbt
)

# Build transcript-to-gene mapping for user annotation
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

# Import and clean user annotation metadata
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

# Handle missing annotation fields and set defaults
# (series of if statements for missing values)

colnames(trann) <- c(
  "gene_id", "gene_biotype", "gene_name",
  "transcript_id", "transcript_biotype", "nearest_ref", "class_code"
)
str_ann <- DataFrame(trann)
trann <- txs_gene_map

# Adjust class codes for transcripts with strand differences
intri <- which(str_ann$class_code == "i")
rnvintri <- unlist(runValue(strand(exs_str[str_ann$transcript_id[intri]])))
rnvintrihost <- unlist(runValue(strand(exonsBy(GTF_annotation, by = "tx", use.names = TRUE)[str_ann$nearest_ref[intri]])))
diffstrand <- names(rnvintri[rnvintri != rnvintrihost])
str_ann[str_ann$transcript_id %in% diffstrand, "class_code"] <- "xi"

# Split user annotation by transcript classification
str_ann_dupl <- str_ann[str_ann[, "class_code"] %in% c("=", "c"), ]
str_ann_new <- str_ann[str_ann[, "class_code"] %in% c("r", "u", "x", "s"), ]
str_ann_isof <- str_ann[str_ann[, "class_code"] %in% c("e", "i", "j", "o", "p", "y", "k"), ]

# Map class codes to descriptive labels
code_annot <- cbind(
  c("=", "c", "k", "r", "u", "x", "s", "e", "i", "j", "o", "p", "y", "m", "n", "xi"),
  c(
    "same", "contained", "extended", "repeat", "intergenic", "antisense", "antisense_spl",
    "pre_mRNA", "intronic", "alt_isof", "exon_overl", "run_on", "enclosing",
    "retain_intr", "retain_partintr", "antisense_intr"
  )
)
code_annot[, 2] <- paste("novel", code_annot[, 2], sep = "_")
str_trann <- data.frame(unique(str_ann[, c("transcript_id", "class_code")]), stringsAsFactors = FALSE)
str_trann[, 2] <- code_annot[match(str_trann[, 2], code_annot[, 1]), 2]

# Set biotype and gene fields for user annotation
str_ann[, "transcript_biotype"] <- "novel"
str_ann$nearest_ref[str_ann$class_code %in% c("r", "u", "x", "s", "xi")] <- NA
mtc <- match(str_ann$nearest_ref, trann$transcript_id)
str_ann[!is.na(mtc), "gene_id"] <- trann$gene_id[mtc[!is.na(mtc)]]
str_ann[!is.na(mtc), "gene_biotype"] <- trann$gene_biotype[mtc[!is.na(mtc)]]
str_ann[!is.na(mtc), "gene_name"] <- trann$gene_name[mtc[!is.na(mtc)]]
str_ann[is.na(mtc), "gene_name"] <- str_ann[is.na(mtc), "gene_id"]
str_ann[is.na(mtc), "gene_biotype"] <- str_ann[is.na(mtc), "transcript_biotype"]

# Merge reference and user annotation based on argument
if (as.logical(args[4]) == FALSE) {
  trann <- rbind(trann, str_ann[!str_ann$class_code %in% c("c", "="), colnames(trann)])
} else if (as.logical(args[4]) == TRUE) {
  trann <- rbind(trann, str_ann[, colnames(trann)])
}

trann$gene_name[is.na(trann$gene_name)] <- trann$gene_id[is.na(trann$gene_name)]

# Combine exons from both annotations and assign metadata
exs_gtf <- unlist(c(exs_str, exs_annot)[trann$transcript_id])
mcols(exs_gtf) <- NULL
exs_gtf$transcript_id <- names(exs_gtf)
exs_gtf$gene_id <- trann[match(exs_gtf$transcript_id, trann$transcript_id), "gene_id"]
exs_gtf$gene_biotype <- trann[match(x = exs_gtf$gene_id, table = trann$gene_id), "gene_biotype"]
exs_gtf$gene_name <- trann[match(x = exs_gtf$gene_id, table = trann$gene_id), "gene_name"]
exs_gtf$transcript_biotype <- trann[match(exs_gtf$transcript_id, trann$transcript_id), "transcript_biotype"]
exs_gtf$type <- "exon"

# Sort and filter exons, ensure strand consistency per gene
all <- sort(exs_gtf)
all$source <- paste0(args[2], "_merged_with_strg")
names(all) <- NULL

all2 <- split(all, all$gene_id)
stra <- strand(all2)
genes_ok <- names(which(elementNROWS(runValue(stra)) == 1))
if (sum(elementNROWS(runValue(stra)) > 1) > 0) {
  stop("Some genes have inconsistent strand info")
}
all_ok <- all2[genes_ok]
all <- sort(unlist(all_ok))
names(all) <- NULL

# Merge STRG and USER gene IDs for transcripts assigned to STRG genes
merged_gtf <- all
user_gtf_txdb <- makeTxDbFromGFF(args[3])
merged_gtf_strg <- merged_gtf[grepl(merged_gtf$transcript_id, pattern = "^R1")]
merged_gtf_user <- merged_gtf[!grepl(merged_gtf$transcript_id, pattern = "^R1")]

# Update gene IDs for user transcripts not assigned to STRG genes
merged_gtf_user$gene_id[!grepl(merged_gtf_user$gene_id, pattern = "^R1")] <-
  paste0("R2.", merged_gtf_user$gene_id[!grepl(merged_gtf_user$gene_id, pattern = "^R1")])

# Map user transcript IDs to user gene IDs
unlist_txs_by_gene <- unlist(transcriptsBy(user_gtf_txdb, by = "gene"))
merged_gtf_user$user_gene_id <-
  names(unlist_txs_by_gene)[match(merged_gtf_user$transcript_id, unlist_txs_by_gene$tx_name)]
merged_gtf_user$user_gene_id <- paste0("R2.", merged_gtf_user$user_gene_id)

# Build mapping between STRG and USER gene IDs for merged genes
df <- as.data.frame(cbind(merged_gtf_user$gene_id, merged_gtf_user$user_gene_id))
df <- df[!duplicated(df), ]
colnames(df) <- c("gene_id", "user_gene_id")
filt_df <- df[grepl(df$gene_id, pattern = "^R1"), ]

# Group gene ID pairs that share values
share_val <- function(row1, row2) {
  any(row1 %in% row2)
}
shares_val <- c(
  FALSE,
  sapply(2:nrow(filt_df), function(i) {
    share_val(filt_df[i, ], filt_df[i - 1, ])
  })
)
filt_df$shares_val <- shares_val
filt_df$group <- cumsum(filt_df$shares_val == FALSE)
split_filt_df <- split(filt_df, filt_df$group)

# Assign new merged gene IDs and write mapping to file
old_gene_ids <- sapply(split_filt_df, function(x) {
  gene_ids <- unique(c(x$gene_id, x$user_gene_id))
  gene_ids <- gene_ids[order(gene_ids)]
  paste(gene_ids, collapse = "-")
})
old_gene_ids <- rep(old_gene_ids, elementNROWS(split_filt_df))
filt_df <- do.call(rbind, split_filt_df)
filt_df$old_gene_ids <- old_gene_ids
filt_df$new_gene_ids <- paste0("R1_R2_merged_", filt_df$group)

gene_id_pairs <- filt_df[, c("old_gene_ids", "new_gene_ids")]
gene_id_pairs <- gene_id_pairs[!duplicated(gene_id_pairs$new_gene_ids), ]
write.table(gene_id_pairs, "gene_ID_pairs.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Update gene IDs in merged user annotation
new_gene_ids <- filt_df$new_gene_ids
names(new_gene_ids) <- filt_df$gene_id
merged_gtf_user_bckp <- merged_gtf_user
merged_gtf_user$gene_id <- new_gene_ids[merged_gtf_user$gene_id]
merged_gtf_user$gene_id[which(is.na(merged_gtf_user$gene_id))] <- merged_gtf_user_bckp$gene_id[which(is.na(merged_gtf_user$gene_id))]
merged_gtf_user$gene_name <- merged_gtf_user$gene_id

# Combine STRG and USER annotations and export merged GTF
merged_gtf <- c(merged_gtf_strg, merged_gtf_user)
merged_gtf <- sort(merged_gtf)

gtf_file <- paste0(args[2], "_merged_with_strg.gtf")
export.gff2(object = merged_gtf, con = gtf_file)
