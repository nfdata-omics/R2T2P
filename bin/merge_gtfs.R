#!/usr/bin/Rscript
suppressPackageStartupMessages(library("RiboseQC"))

# takes standard and de-novo GTFs
args <- commandArgs(trailingOnly = TRUE)

GTF_annotation <- get(load(args[1]))
annotation_ok <- GTF_annotation
annotation_str <- makeTxDbFromGFF(file = args[2])

exs_annot <- GTF_annotation$exons_txs
exs_str <- exonsBy(annotation_str, by = "tx", use.names = T)
genes_annot <- GTF_annotation$genes
genes_str <- GenomicFeatures::genes(annotation_str)

txs_gene_map <- GTF_annotation$trann

transcripts_db_str <- transcripts(annotation_str)
txs_gene_map_str <- mapIds(transcripts_db_str$tx_name, x = annotation_str, column = "GENEID", keytype = "TXNAME")
txs_gene_map_str <- data.frame(gene_id = txs_gene_map_str, transcript_id = names(txs_gene_map_str), stringsAsFactors = F)
rownames(txs_gene_map_str) <- NULL

# take away the exact same txs and adjust the gene ids of the new ones,
# based on overlaps with annotated genes
trann <- unique(mcols(import.gff2(args[2], colnames = c(
  "gene_id",
  "gene_biotype", "gene_type", "gene_name", "gene_symbol",
  "transcript_id", "transcript_biotype", "transcript_type",
  "cmp_ref", "class_code"
))))
trann <- trann[!is.na(trann$class_code), ]
trann <- data.frame(unique(trann), stringsAsFactors = F)
if (sum(!is.na(trann$transcript_biotype)) == 0 & sum(!is.na(trann$transcript_type)) == 0) {
  trann$transcript_biotype <- "no_type"
}
if (sum(!is.na(trann$transcript_biotype)) == 0) {
  trann$transcript_biotype <- NULL
}
if (sum(!is.na(trann$transcript_type)) == 0) {
  trann$transcript_type <- NULL
}
if (sum(!is.na(trann$gene_biotype)) == 0 & sum(!is.na(trann$gene_type)) == 0) {
  trann$gene_type <- "no_type"
}
if (sum(!is.na(trann$gene_name)) == 0 & sum(!is.na(trann$gene_symbol)) == 0) {
  trann$gene_name <- "no_name"
}
if (sum(!is.na(trann$gene_biotype)) == 0) {
  trann$gene_biotype <- NULL
}
if (sum(!is.na(trann$gene_type)) == 0) {
  trann$gene_type <- NULL
}
if (sum(!is.na(trann$gene_name)) == 0) {
  trann$gene_name <- NULL
}
if (sum(!is.na(trann$gene_symbol)) == 0) {
  trann$gene_symbol <- NULL
}
colnames(trann) <- c(
  "gene_id", "gene_biotype", "gene_name",
  "transcript_id", "transcript_biotype", "nearest_ref", "class_code"
)
str_ann <- DataFrame(trann)
trann <- GTF_annotation$trann
intri <- which(str_ann$class_code == "i")

rnvintri <- unlist(runValue(strand(exs_str[str_ann$transcript_id[intri]])))
rnvintrihost <- unlist(runValue(strand(GTF_annotation$exons_txs[str_ann$nearest_ref[intri]])))
diffstrand <- names(rnvintri[rnvintri != rnvintrihost])
str_ann[str_ann$transcript_id %in% diffstrand, "class_code"] <- "xi"

str_ann_dupl <- str_ann[str_ann[, "class_code"] %in% c("=", "c"), ]
str_ann_new <- str_ann[str_ann[, "class_code"] %in% c("r", "u", "x", "s"), ]
str_ann_isof <- str_ann[str_ann[, "class_code"] %in% c("e", "i", "j", "o", "p", "y", "k"), ]
code_annot <- cbind(c("=", "c", "k", "r", "u", "x", "s", "e", "i", "j", "o", "p", "y", "m", "n", "xi"), c("same", "contained", "extended", "repeat", "intergenic", "antisense", "antisense_spl", "pre_mRNA", "intronic", "alt_isof", "exon_overl", "run_on", "enclosing", "retain_intr", "retain_partintr", "antisense_intr"))
code_annot[, 2] <- paste("novel", code_annot[, 2], sep = "_")
str_trann <- data.frame(unique(str_ann[, c("transcript_id", "class_code")]), stringsAsFactors = F)
str_trann[, 2] <- code_annot[match(str_trann[, 2], code_annot[, 1]), 2]

str_ann[, "transcript_biotype"] <- code_annot[match(str_ann$class_code, code_annot[, 1]), 2]
str_ann$nearest_ref[str_ann$class_code %in% c("r", "u", "x", "s", "xi")] <- NA
mtc <- match(str_ann$nearest_ref, trann$transcript_id)
str_ann[!is.na(mtc), "gene_id"] <- trann$gene_id[mtc[!is.na(mtc)]]
str_ann[!is.na(mtc), "gene_biotype"] <- trann$gene_biotype[mtc[!is.na(mtc)]]
str_ann[!is.na(mtc), "gene_name"] <- trann$gene_name[mtc[!is.na(mtc)]]
str_ann[is.na(mtc), "gene_name"] <- str_ann[is.na(mtc), "gene_id"]
str_ann[is.na(mtc), "gene_biotype"] <- str_ann[is.na(mtc), "transcript_biotype"]

if (as.logical(args[3]) == FALSE) {
  trann <- rbind(trann, str_ann[!str_ann$class_code %in% c("c", "="), colnames(trann)])
} else if (as.logical(args[3]) == TRUE) {
  trann <- rbind(trann, str_ann[, colnames(trann)])
}

trann$gene_name[is.na(trann$gene_name)] <- trann$gene_id[is.na(trann$gene_name)]

exs_gtf <- unlist(c(exs_str, exs_annot)[trann$transcript_id])
mcols(exs_gtf) <- NULL
exs_gtf$transcript_id <- names(exs_gtf)
exs_gtf$gene_id <- trann[match(exs_gtf$transcript_id, trann$transcript_id), "gene_id"]
exs_gtf$gene_biotype <- trann[match(x = exs_gtf$gene_id, table = trann$gene_id), "gene_biotype"]
exs_gtf$gene_name <- trann[match(x = exs_gtf$gene_id, table = trann$gene_id), "gene_name"]
exs_gtf$transcript_biotype <- trann[match(exs_gtf$transcript_id, trann$transcript_id), "transcript_biotype"]
exs_gtf$type <- "exon"

cds_gtf <- unlist(GTF_annotation$cds_txs)
mcols(cds_gtf) <- NULL
cds_gtf$transcript_id <- names(cds_gtf)
cds_gtf$gene_id <- trann[match(cds_gtf$transcript_id, trann$transcript_id), "gene_id"]
cds_gtf$gene_biotype <- trann[match(x = cds_gtf$gene_id, table = trann$gene_id), "gene_biotype"]
cds_gtf$gene_name <- trann[match(x = cds_gtf$gene_id, table = trann$gene_id), "gene_name"]
cds_gtf$transcript_biotype <- trann[match(cds_gtf$transcript_id, trann$transcript_id), "transcript_biotype"]
cds_gtf$type <- "CDS"

all <- sort(c(exs_gtf, cds_gtf))
all$`source` <- paste0(args[2], "_merged_with_ref")
names(all) <- NULL

all2 <- split(all, all$gene_id)
stra <- strand(all2)
genes_ok <- names(which(elementNROWS(runValue(stra)) == 1))
if (sum(elementNROWS(runValue(stra)) > 1) > 0) {
  stop("Some genes have incosistent strand info")
}
all_ok <- all2[genes_ok]

all <- sort(unlist(all_ok))
names(all) <- NULL

gtf_file <- paste0(args[2], "_stringtie.gtf")

export.gff2(object = all, con = gtf_file)
gtf_lines <- readLines(gtf_file, warn = FALSE)
writeLines(gtf_lines[!grepl("^##date ", gtf_lines)], gtf_file)

seqinfotwob <- GTF_annotation$seqinfo
annotation <- makeTxDbFromGFF(file = gtf_file, format = "gtf", chrominfo = seqinfotwob)

genes <- GenomicFeatures::genes(annotation)
exons_ge <- exonsBy(annotation, by = "gene")
exons_ge <- reduce(exons_ge)

cds_gen <- cdsBy(annotation, "gene")
cds_ge <- reduce(cds_gen)

# define regions not overlapping CDS ( or exons when defining introns)
threeutrs <- reduce(GenomicRanges::setdiff(unlist(threeUTRsByTranscript(annotation)), unlist(cds_ge), ignore.strand = FALSE))
fiveutrs <- reduce(GenomicRanges::setdiff(unlist(fiveUTRsByTranscript(annotation)), unlist(cds_ge), ignore.strand = FALSE))
introns <- reduce(GenomicRanges::setdiff(unlist(intronsByTranscript(annotation)), unlist(exons_ge), ignore.strand = FALSE))
nc_exons <- reduce(GenomicRanges::setdiff(unlist(exons_ge), reduce(c(unlist(cds_ge), fiveutrs, threeutrs)), ignore.strand = FALSE))

# assign gene ids (mutiple when overlapping multiple genes)
ov <- findOverlaps(threeutrs, genes)
ov <- split(subjectHits(ov), queryHits(ov))
threeutrs$gene_id <- CharacterList(lapply(ov, FUN = function(x) {
  names(genes)[x]
}))
ov <- findOverlaps(fiveutrs, genes)
ov <- split(subjectHits(ov), queryHits(ov))
fiveutrs$gene_id <- CharacterList(lapply(ov, FUN = function(x) {
  names(genes)[x]
}))
ov <- findOverlaps(introns, genes)
ov <- split(subjectHits(ov), queryHits(ov))
introns$gene_id <- CharacterList(lapply(ov, FUN = function(x) {
  names(genes)[x]
}))
ov <- findOverlaps(nc_exons, genes)
ov <- split(subjectHits(ov), queryHits(ov))
nc_exons$gene_id <- CharacterList(lapply(ov, FUN = function(x) {
  names(genes)[x]
}))

intergenicRegions <- genes
strand(intergenicRegions) <- "*"
intergenicRegions <- gaps(reduce(intergenicRegions))
intergenicRegions <- intergenicRegions[strand(intergenicRegions) == "*"]

cds_tx <- cdsBy(annotation, "tx", use.names = T)
txs_gene <- transcriptsBy(annotation, by = "gene")
genes_red <- reduce(sort(GenomicFeatures::genes(annotation)))

exons_tx <- exonsBy(annotation, "tx", use.names = T)

transcripts_db <- transcripts(annotation)
intron_names_tx <- intronsByTranscript(annotation, use.names = T)

# define exonic bins, including regions overlapping multiple genes
nsns <- exonicParts(annotation, linked.to.single.gene.only = F)

# define tx_coordinates of ORF boundaries
exsss_cds <- exons_tx[names(cds_tx)]
chunks <- seq(1, length(cds_tx), by = 20000)
if (chunks[length(chunks)] < length(cds_tx)) {
  chunks <- c(chunks, length(cds_tx))
}
mapp <- GRangesList()
for (i in 1:(length(chunks) - 1)) {
  if (i != (length(chunks) - 1)) {
    mapp <- suppressWarnings(c(mapp, pmapToTranscripts(cds_tx[chunks[i]:(chunks[i + 1] - 1)], transcripts = exsss_cds[chunks[i]:(chunks[i + 1] - 1)])))
  }
  if (i == (length(chunks) - 1)) {
    mapp <- suppressWarnings(c(mapp, pmapToTranscripts(cds_tx[chunks[i]:(chunks[i + 1])], transcripts = exsss_cds[chunks[i]:(chunks[i + 1])])))
  }
}
cds_txscoords <- unlist(mapp)

# extract biotypes and ids
trann <- unique(mcols(import.gff2(gtf_file, colnames = c("gene_id", "gene_biotype", "gene_type", "gene_name", "gene_symbol", "transcript_id", "transcript_biotype", "transcript_type"))))
trann <- trann[!is.na(trann$transcript_id), ]
trann <- data.frame(unique(trann), stringsAsFactors = F)

if (sum(!is.na(trann$transcript_biotype)) == 0 & sum(!is.na(trann$transcript_type)) == 0) {
  trann$transcript_biotype <- "no_type"
}
if (sum(!is.na(trann$transcript_biotype)) == 0) {
  trann$transcript_biotype <- NULL
}
if (sum(!is.na(trann$transcript_type)) == 0) {
  trann$transcript_type <- NULL
}

if (sum(!is.na(trann$gene_biotype)) == 0 & sum(!is.na(trann$gene_type)) == 0) {
  trann$gene_type <- "no_type"
}
if (sum(!is.na(trann$gene_name)) == 0 & sum(!is.na(trann$gene_symbol)) == 0) {
  trann$gene_name <- "no_name"
}
if (sum(!is.na(trann$gene_biotype)) == 0) {
  trann$gene_biotype <- NULL
}
if (sum(!is.na(trann$gene_type)) == 0) {
  trann$gene_type <- NULL
}
if (sum(!is.na(trann$gene_name)) == 0) {
  trann$gene_name <- NULL
}
if (sum(!is.na(trann$gene_symbol)) == 0) {
  trann$gene_symbol <- NULL
}
colnames(trann) <- c("gene_id", "gene_biotype", "gene_name", "transcript_id", "transcript_biotype")

trann <- DataFrame(trann)

# introns and transcript_ids/gene_ids
unq_intr <- sort(unique(unlist(intron_names_tx)))
names(unq_intr) <- NULL
all_intr <- unlist(intron_names_tx)

ov <- findOverlaps(unq_intr, all_intr, type = "equal")
ov <- split(subjectHits(ov), queryHits(ov))
a_nam <- CharacterList(lapply(ov, FUN = function(x) {
  unique(names(all_intr)[x])
}))

unq_intr$type <- "J"
unq_intr$tx_name <- a_nam

mat_genes <- match(unq_intr$tx_name, trann$transcript_id)
g <- unlist(apply(cbind(1:length(mat_genes), Y = elementNROWS(mat_genes)), FUN = function(x) rep(x[1], x[2]), MARGIN = 1))
g2 <- split(trann[unlist(mat_genes), "gene_id"], g)
unq_intr$gene_id <- CharacterList(lapply(g2, unique))

# filter ncRNA and ncIsof regions
ncrnas <- nc_exons[!nc_exons %over% genes[trann$gene_id[trann$gene_biotype == "protein_coding"]]]
ncisof <- nc_exons[nc_exons %over% genes[trann$gene_id[trann$gene_biotype == "protein_coding"]]]

cds_txscoords$gene_id <- trann$gene_id[match(as.vector(seqnames(cds_txscoords)), trann$transcript_id)]
cds_cc <- cds_txscoords
strand(cds_cc) <- "*"
sta_cc <- resize(cds_cc, width = 1, "start")
sta_cc <- unlist(pmapFromTranscripts(sta_cc, exons_tx[seqnames(sta_cc)], ignore.strand = F))
sta_cc$gene_id <- trann$gene_id[match(names(sta_cc), trann$transcript_id)]
sta_cc <- sta_cc[sta_cc$hit]
strand(sta_cc) <- structure(as.vector(strand(transcripts_db)), names = transcripts_db$tx_name)[names(sta_cc)]
sta_cc$type <- "start_codon"
mcols(sta_cc) <- mcols(sta_cc)[, c("exon_rank", "type", "gene_id")]

sto_cc <- resize(cds_cc, width = 1, "end")
# stop codon is the 1st nt, e.g. U of the UAA
sto_cc <- shift(sto_cc, -2)
stop_inannot <- GTF_annotation$stop_in_gtf
if (is.na(stop_inannot)) {
  sto_cc <- resize(trim(shift(sto_cc, 3)), width = 1, fix = "end")
}

sto_cc <- unlist(pmapFromTranscripts(sto_cc, exons_tx[seqnames(sto_cc)], ignore.strand = F))
sto_cc <- sto_cc[sto_cc$hit]
sto_cc$gene_id <- trann$gene_id[match(names(sto_cc), trann$transcript_id)]
strand(sto_cc) <- structure(as.vector(strand(transcripts_db)), names = transcripts_db$tx_name)[names(sto_cc)]
sto_cc$type <- "stop_codon"
mcols(sto_cc) <- mcols(sto_cc)[, c("exon_rank", "type", "gene_id")]

# define most common, most upstream/downstream
start_stop_cc <- sort(c(sta_cc, sto_cc))
start_stop_cc$transcript_id <- names(start_stop_cc)
start_stop_cc$most_up_downstream <- FALSE
start_stop_cc$most_frequent <- FALSE

df <- cbind.DataFrame(start(start_stop_cc), start_stop_cc$type, start_stop_cc$gene_id)
colnames(df) <- c("start_pos", "type", "gene_id")
upst <- by(df$start_pos, INDICES = df$gene_id, function(x) {
  x == min(x) | x == max(x)
})
start_stop_cc$most_up_downstream <- unlist(upst[unique(df$gene_id)])

mostfr <- by(df[, c("start_pos", "type")], INDICES = df$gene_id, function(x) {
  mfreq <- table(x)
  x$start_pos %in% as.numeric(names(which(mfreq[, 1] == max(mfreq[, 1])))) | x$start_pos %in% as.numeric(names(which(mfreq[, 2] == max(mfreq[, 2]))))
})

start_stop_cc$most_frequent <- unlist(mostfr[unique(df$gene_id)])

names(start_stop_cc) <- NULL

# define transcripts as containing frequent start/stop codons or most upstream ones, in relation with 5'UTR length
mostupstr_tx <- sum(LogicalList(split(start_stop_cc$most_up_downstream, start_stop_cc$transcript_id)))[as.character(seqnames(cds_txscoords))]
cds_txscoords$upstr_stasto <- mostupstr_tx
mostfreq_tx <- sum(LogicalList(split(start_stop_cc$most_frequent, start_stop_cc$transcript_id)))[as.character(seqnames(cds_txscoords))]
cds_txscoords$mostfreq_stasto <- mostfreq_tx
cds_txscoords$lentx <- sum(width(exons_tx[as.character(seqnames(cds_txscoords))]))
df <- cbind.DataFrame(as.character(seqnames(cds_txscoords)), width(cds_txscoords), start(cds_txscoords), cds_txscoords$mostfreq_stasto, cds_txscoords$gene_id)
colnames(df) <- c("txid", "cdslen", "utr5len", "var", "gene_id")
repres_freq <- by(df[, c("txid", "cdslen", "utr5len", "var")], df$gene_id, function(x) {
  x <- x[order(x$var, x$utr5len, x$cdslen, decreasing = T), ]
  x <- x[x$var == max(x$var), ]
  ok <- x$txid[which(x$cdslen == max(x$cdslen) & x$utr5len == max(x$utr5len) & x$var == max(x$var))][1]
  if (length(ok) == 0 | is.na(ok[1])) {
    ok <- x$txid[1]
  }
  ok
})

df <- cbind.DataFrame(as.character(seqnames(cds_txscoords)), width(cds_txscoords), start(cds_txscoords), cds_txscoords$upstr_stasto, cds_txscoords$gene_id)
colnames(df) <- c("txid", "cdslen", "utr5len", "var", "gene_id")
repres_upstr <- by(df[, c("txid", "cdslen", "utr5len", "var")], df$gene_id, function(x) {
  x <- x[order(x$var, x$utr5len, x$utr5len, decreasing = T), ]
  x <- x[x$var == max(x$var), ]
  ok <- x$txid[which(x$cdslen == max(x$cdslen) & x$utr5len == max(x$utr5len) & x$var == max(x$var))][1]
  if (length(ok) == 0 | is.na(ok[1])) {
    ok <- x$txid[1]
  }
  ok
})
df <- cbind.DataFrame(as.character(seqnames(cds_txscoords)), width(cds_txscoords), start(cds_txscoords), cds_txscoords$upstr_stasto, cds_txscoords$gene_id)
colnames(df) <- c("txid", "cdslen", "utr5len", "var", "gene_id")
repres_len5 <- by(df[, c("txid", "cdslen", "utr5len", "var")], df$gene_id, function(x) {
  x <- x[order(x$utr5len, x$var, x$cdslen, decreasing = T), ]
  ok <- x$txid[which(x$utr5len == max(x$utr5len) & x$var == max(x$var))][1]
  if (length(ok) == 0 | is.na(ok[1])) {
    ok <- x$txid[1]
  }
  ok
})

cds_txscoords$reprentative_mostcommon <- as.character(seqnames(cds_txscoords)) %in% unlist(repres_freq)
cds_txscoords$reprentative_boundaries <- as.character(seqnames(cds_txscoords)) %in% unlist(repres_upstr)
cds_txscoords$reprentative_5len <- as.character(seqnames(cds_txscoords)) %in% unlist(repres_len5)
unq_stst <- start_stop_cc
mcols(unq_stst) <- NULL
unq_stst <- sort(unique(unq_stst))
ov <- findOverlaps(unq_stst, start_stop_cc, type = "equal")
ov <- split(subjectHits(ov), queryHits(ov))
unq_stst$type <- CharacterList(lapply(ov, FUN = function(x) {
  unique(start_stop_cc$type[x])
}))
unq_stst$transcript_id <- CharacterList(lapply(ov, FUN = function(x) {
  start_stop_cc$transcript_id[x]
}))
unq_stst$gene_id <- CharacterList(lapply(ov, FUN = function(x) {
  unique(start_stop_cc$gene_id[x])
}))

unq_stst$reprentative_mostcommon <- sum(!is.na(match(unq_stst$transcript_id, unlist(as(repres_freq, "CharacterList"))))) > 0
unq_stst$reprentative_boundaries <- sum(!is.na(match(unq_stst$transcript_id, unlist(as(repres_upstr, "CharacterList"))))) > 0
unq_stst$reprentative_5len <- sum(!is.na(match(unq_stst$transcript_id, unlist(as(repres_len5, "CharacterList"))))) > 0

# put in a list
GTF_annotation <- list(transcripts_db, txs_gene, GTF_annotation$seqinfo, unq_stst, cds_tx, intron_names_tx, cds_gen, exons_tx, nsns, unq_intr, genes, threeutrs, fiveutrs, ncisof, ncrnas, introns, intergenicRegions, trann, cds_txscoords, GTF_annotation$genetic_codes, GTF_annotation$genome_package, GTF_annotation$stop_in_gtf)
names(GTF_annotation) <- c("txs", "txs_gene", "seqinfo", "start_stop_codons", "cds_txs", "introns_txs", "cds_genes", "exons_txs", "exons_bins", "junctions", "genes", "threeutrs", "fiveutrs", "ncIsof", "ncRNAs", "introns", "intergenicRegions", "trann", "cds_txs_coords", "genetic_codes", "genome_package", "stop_in_gtf")

txs_all <- unique(GTF_annotation$trann$transcript_id)
txs_exss <- unique(names(GTF_annotation$exons_txs))

txs_notok <- txs_all[!txs_all %in% txs_exss]
if (length(txs_notok) > 0) {
  set.seed(666)
  cat(paste("Warning: ", length(txs_notok), " txs with incorrect/unspecified exon boundaries - e.g. trans-splicing events, examples: ",
    paste(txs_notok[sample(1:length(txs_notok), size = min(3, length(txs_notok)), replace = F)], collapse = ", "), " - ", date(), "\n",
    sep = ""
  ))
}

rannot_file <- paste0(args[2], "_stringtie_Rannot")
save(GTF_annotation, file = rannot_file)
