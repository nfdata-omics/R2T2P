# Attaching R packages
suppressMessages(library(cowplot))
suppressMessages(library(zoo))
suppressMessages(library(reshape2))
suppressMessages(library(limma))
suppressMessages(library(limpa))
suppressMessages(library(rtracklayer))
suppressMessages(library(ORFquant))

# Command-line arguments.
cli_args <- commandArgs(trailingOnly = TRUE)

print_usage <- function() {
  cat("Usage: Rscript diff_proteomics.R [options]\n\n",
      "  --data-type <DDA_LFQ|DDA_TMT|DIA>\n",
      "  --gtf-rannot <RData file>\n",
      "  --orfquant-results <RData file>\n",
      "  --search-fold <FragPipe results directory>\n",
      "  --manifest <file>\n",
      "  --baseline <condition label>\n",
      "  --annotation-files <comma-separated files>\n",
      "  --help\n",
      sep = "")
}

get_cli_arg <- function(option, default = NULL, required = FALSE) {
  equals_pattern <- paste0("^", option, "=")
  equals_match <- grep(equals_pattern, cli_args, value = TRUE)
  if (length(equals_match) > 0) return(sub(equals_pattern, "", equals_match[[1]]))

  position <- match(option, cli_args)
  if (!is.na(position)) {
    if (position == length(cli_args) || startsWith(cli_args[[position + 1]], "--")) {
      stop("Missing value for ", option, call. = FALSE)
    }
    return(cli_args[[position + 1]])
  }
  if (required) stop("Missing required argument ", option, call. = FALSE)
  default
}

if ("--help" %in% cli_args) {
  print_usage()
  quit(status = 0)
}

if (length(cli_args) > 0) {
  data_type <- get_cli_arg("--data-type", required = TRUE)
  if (!data_type %in% c("DDA_LFQ", "DDA_TMT", "DIA")) {
    stop("--data-type must be one of: DDA_LFQ, DDA_TMT, DIA", call. = FALSE)
  }
  gtf_rannot_path <- get_cli_arg("--gtf-rannot", required = TRUE)
  orfquant_results_path <- get_cli_arg("--orfquant-results", required = TRUE)
  search_fold <- get_cli_arg("--search-fold", required = TRUE)
  combined_peptide_path <- file.path(search_fold, "combined_peptide.tsv")
  abundance_peptide_path <- file.path(search_fold, "tmt-report", "abundance_peptide_None.tsv")
  dia_abundance_path <- file.path(search_fold, "dia-quant-output", "abundance_peptide_MS2quant_None.tsv")
  manifest_path <- get_cli_arg("--manifest")
  baseline <- get_cli_arg("--baseline", required = TRUE)

  annotation_files <- get_cli_arg("--annotation-files")
  annotation_fls_paths <- if (is.null(annotation_files)) character() else strsplit(annotation_files, ",", fixed = TRUE)[[1]]

  GTF_annotation <- get(load(gtf_rannot_path))
  load(orfquant_results_path)

  if (data_type %in% c("DDA_LFQ", "DDA_TMT") && !file.exists(combined_peptide_path)) {
    stop("Missing combined_peptide.tsv: ", combined_peptide_path, call. = FALSE)
  }
  if (data_type == "DIA" && !file.exists(dia_abundance_path)) {
    stop("Missing DIA abundance file: ", dia_abundance_path, call. = FALSE)
  }
  if (data_type == "DDA_TMT" && !file.exists(abundance_peptide_path)) {
    stop("Missing TMT abundance file: ", abundance_peptide_path, call. = FALSE)
  }
  if (data_type == "DDA_TMT" && length(annotation_fls_paths) == 0) {
    stop("--annotation-files is required when --data-type DDA_TMT", call. = FALSE)
  }
  if (data_type != "DDA_TMT" && is.null(manifest_path)) {
    stop("--manifest is required when --data-type is DDA_LFQ or DIA", call. = FALSE)
  }
}

# The "combined_peptide.tsv" file (or the "abundance_peptide_MS2quant_None.tsv" file for DIA)
if (data_type == "DIA") {
  combined_peptide <- read.table(dia_abundance_path, sep="\t", header=T)
} else {
  combined_peptide <- read.table(combined_peptide_path, sep="\t", header=T)
}
combined_peptide_prots <- CharacterList(strsplit(paste(combined_peptide$Protein, combined_peptide$Mapped.Proteins, sep=", "), split=", "))
names(combined_peptide_prots) <- combined_peptide$Peptide.Sequence
keep <- !sapply(combined_peptide_prots, function(x){any(grepl(pattern="^(rev|contam)",x))})
peps_to_remove <- names(combined_peptide_prots)[!keep]
combined_peptide_prots <- combined_peptide_prots[keep]
combined_peptide <- combined_peptide[combined_peptide$Peptide.Sequence %in% names(combined_peptide_prots),]
rownames(combined_peptide) <- combined_peptide$Peptide.Sequence

if (data_type == "DDA_TMT"){

  # The "abundance_peptide_[normalization].tsv" file
  abundance_peptide <- read.table(abundance_peptide_path, sep="\t", header=T)
  rownames(abundance_peptide) <- abundance_peptide$Peptide
  abundance_peptide <- abundance_peptide[rownames(abundance_peptide) %in% rownames(combined_peptide),]
  combined_peptide <- combined_peptide[rownames(combined_peptide) %in% rownames(abundance_peptide),]
  abundance_peptide <- abundance_peptide[rownames(combined_peptide),]

  # Mapping peptides to proteins
  abundance_peptide_prots <- CharacterList(strsplit(paste(abundance_peptide$Protein, abundance_peptide$Mapped.Proteins, sep=", "), split=", "))
  names(abundance_peptide_prots) <- abundance_peptide$Peptide
  peps_to_remove <- unique(c(peps_to_remove, names(abundance_peptide_prots)[sapply(abundance_peptide_prots, function(x) {any(grepl(pattern="^(rev|contam)",x))})]))
  combined_peptide <- combined_peptide[!(abundance_peptide$Peptide %in% peps_to_remove),]
  abundance_peptide <- abundance_peptide[!(abundance_peptide$Peptide %in% peps_to_remove),]
  peps_mapinfo <- DataFrame(abundance_peptide[,c("Protein","ProteinID","Entry.Name","Protein.Description","Mapped.Proteins","Gene","Mapped.Genes","Start","MaxPepProb")])

} else {

  # Mapping peptides to proteins
  peps_mapinfo <- DataFrame(combined_peptide[,c("Protein","Protein.ID","Entry.Name","Protein.Description","Mapped.Proteins","Gene","Mapped.Genes","Start")])

}

peps_mapinfo$All.Mapped.Proteins <- CharacterList(strsplit(paste(peps_mapinfo$Protein, peps_mapinfo$Mapped.Proteins, sep=", "), split=", "))
peptide_to_proteins <- peps_mapinfo$All.Mapped.Proteins
names(peptide_to_proteins) <- rownames(peps_mapinfo)

# Expanding and adding columns (gene_id, adj_unique_features_reads, ORF_category_Gen, ORF_category_Tx_compatible, and compatible_biotype)
exp_peptide_to_proteins <- S4Vectors::expand(DataFrame(peptide_to_proteins), colnames=names(DataFrame(peptide_to_proteins)))
exp_peptide_to_proteins$transcript_id <- sapply(strsplit(gsub(exp_peptide_to_proteins$peptide_to_proteins, pattern="^R2T2P_", replacement=""), split="_"), "[[", 1)
exp_peptide_to_proteins$gene_id <- GTF_annotation$trann$gene_id[match(exp_peptide_to_proteins$transcript_id, GTF_annotation$trann$transcript_id)]
exp_peptide_to_proteins$transcript_id <- NULL
adj_unique_features_reads <- sum(ORFquant_results$ORFs_tx$adj_unique_features_reads)
ORF_category_Gen <- ORFquant_results$ORFs_tx$ORF_category_Gen
ORF_category_Tx_compatible <- ORFquant_results$ORFs_tx$ORF_category_Tx_compatible
orf_ids <- ORFquant_results$ORFs_tx$ORF_id_tr
ORFs_tx<-ORFquant_results$ORFs_tx
ORFs_tx$tx_novel<-sum(!grepl(ORFs_tx$compatible_with,pattern="^R1|^R2"))==0
ORFs_tx$tx_novel_R1R2<-"none"
ORFs_tx$tx_novel_R1R2[ORFs_tx$tx_novel]<-"R1R2"
ORFs_tx$tx_novel_R1R2[sum(!grepl(ORFs_tx$compatible_with,pattern="^R1"))==0]<-"R1"
ORFs_tx$tx_novel_R1R2[sum(!grepl(ORFs_tx$compatible_with,pattern="^R2"))==0]<-"R2"
ORFs_tx$compatible_biotype[ORFs_tx$tx_novel]<-paste("novel_tx",ORFs_tx$tx_novel_R1R2[ORFs_tx$tx_novel],sep="_")
compatible_biotype<-ORFs_tx$compatible_biotype
names(compatible_biotype) <- orf_ids
names(adj_unique_features_reads) <- orf_ids
names(ORF_category_Gen) <- orf_ids
names(ORF_category_Tx_compatible) <- orf_ids
exp_peptide_to_proteins_orfq <- exp_peptide_to_proteins[grepl(exp_peptide_to_proteins$peptide_to_proteins, pattern="^R2T2P_"),]
orfq_orf_ids <- gsub(exp_peptide_to_proteins_orfq$peptide_to_proteins, pattern="^R2T2P_", replacement="")
exp_peptide_to_proteins_orfq$adj_unique_features_reads <- adj_unique_features_reads[orfq_orf_ids]
exp_peptide_to_proteins_orfq$ORF_category_Gen <- ORF_category_Gen[orfq_orf_ids]
exp_peptide_to_proteins_orfq$ORF_category_Tx_compatible <- ORF_category_Tx_compatible[orfq_orf_ids]
exp_peptide_to_proteins_orfq$compatible_biotype <- compatible_biotype[orfq_orf_ids]
if (search_fold != "custom_database_search"){
  exp_peptide_to_proteins_genc <- exp_peptide_to_proteins[!grepl(exp_peptide_to_proteins$peptide_to_proteins, pattern="^R2T2P_"),]
  exp_peptide_to_proteins_genc$adj_unique_features_reads <- NA
  exp_peptide_to_proteins_genc$ORF_category_Gen <- NA
  exp_peptide_to_proteins_genc$ORF_category_Tx_compatible <- NA
  exp_peptide_to_proteins_genc$compatible_biotype <- NA
  exp_peptide_to_proteins <- rbind(exp_peptide_to_proteins_genc, exp_peptide_to_proteins_orfq)
} else {
  exp_peptide_to_proteins <- exp_peptide_to_proteins_orfq
}

group_by_peptide <- function(x, numeric = FALSE){
  x <- split(x, rownames(exp_peptide_to_proteins))
  if (numeric) unname(NumericList(x)) else unname(CharacterList(x))
}

# Peptide annotation
peptide_annotation <- DataFrame(Mapped.Proteins = group_by_peptide(exp_peptide_to_proteins$peptide_to_proteins),
                                Mapped.Genes = group_by_peptide(exp_peptide_to_proteins$gene_id),
                                adj_unique_features_reads = group_by_peptide(exp_peptide_to_proteins$adj_unique_features_reads),
                                ORF_category_Gen = group_by_peptide(exp_peptide_to_proteins$ORF_category_Gen),
                                ORF_category_Tx_compatible = group_by_peptide(exp_peptide_to_proteins$ORF_category_Tx_compatible),
                                compatible_biotype = group_by_peptide(exp_peptide_to_proteins$compatible_biotype))
rownames(peptide_annotation) <- names(CharacterList(split(exp_peptide_to_proteins$peptide_to_proteins, rownames(exp_peptide_to_proteins))))
peptide_annotation$MainProtein <- peps_mapinfo$Protein[match(rownames(peptide_annotation), rownames(peps_mapinfo))]
peptide_annotation$Start_on_MainProtein <- peps_mapinfo$Start[match(rownames(peptide_annotation), rownames(peps_mapinfo))]

if (data_type == "DDA_TMT"){

  # Table for DE
  annots_df <- do.call("rbind", lapply(annotation_fls_paths, function(x){read.table(x, sep=" ")}))
  colnames(annots_df) <- c("channel", "sample_name")
  annots_df <- annots_df[!is.na(annots_df$sample_name) & annots_df$sample_name != "NA" & !grepl(pattern="^Reference", annots_df$sample_name),]
  de_table <- annots_df
  de_table$condition <- sapply(strsplit(de_table$sample_name, "_(?=[^_]+$)", perl = TRUE), `[`, 1)
  de_table$replicate <- sapply(strsplit(de_table$sample_name, "_(?=[^_]+$)", perl = TRUE), `[`, 2)
  de_table$baseline <- (de_table$condition == baseline)
  de_table <- de_table[,c("condition", "replicate", "sample_name", "baseline")]
  samples <- de_table$sample_name

} else {

  # Table for DE
  manifest <- read.table(manifest_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  colnames(manifest) <- c("Path", "Experiment", "Bioreplicate", "Data type")
  # -- Extended table [one row per file]
  de_table <- manifest[,c("Path", "Experiment", "Bioreplicate")]
  colnames(de_table) <- c("file_path", "condition", "replicate")
  de_table$sample_name <- paste(de_table$condition, de_table$replicate, sep="_")
  de_table$baseline <- (de_table$condition == baseline)
  # -- Simplified table [one row per sample, i.e., condition + replicate]
  simpl_de_table <- de_table[!duplicated(de_table$sample_name),]
  simpl_de_table$file_path <- NULL

  if (data_type == "DDA_LFQ") {
    if (any(grepl(pattern="MaxLFQ.Intensity$", colnames(combined_peptide)))) {
      abundance_peptide <- combined_peptide[,grepl(pattern="MaxLFQ.Intensity$", colnames(combined_peptide))]
      samples <- paste0(simpl_de_table$sample_name, ".MaxLFQ.Intensity")
    } else {
      abundance_peptide <- combined_peptide[,grepl(pattern="Intensity$", colnames(combined_peptide))]
      samples <- paste0(simpl_de_table$sample_name, ".Intensity")
    }
  } else if (data_type == "DIA") {
    abundance_peptide <- combined_peptide[,16:ncol(combined_peptide)]
    samples <- sub("\\.[^.]+$", "", basename(de_table$file_path))
  }

}

# Peptide intensities
sampok <- samples %in% colnames(abundance_peptide)
if(sum(!sampok)>0){
  warning(paste("Samples not found in intensity table:", paste(samples[!sampok], collapse="; ")))
}
sampok <- samples[sampok]
peptide_intensities <- abundance_peptide[,colnames(abundance_peptide) %in% sampok]
peptide_intensities <- peptide_intensities[,sampok]

# Saving peptide annotation as an RData file
save(peptide_annotation, file = paste0(search_fold, "_peptide_annotation.RData"))
# Writing peptide intensities to a tsv file
write.table(peptide_intensities, file = paste0(search_fold, "_peptide_intensities.tsv"), sep = "\t", quote = F, row.names = T, col.names = T)

# Metadata
df_meta <- peptide_annotation
colnames(df_meta)[colnames(df_meta)=="Mapped.Genes"]<-"Mapped_Genes"
colnames(df_meta)[colnames(df_meta)=="Mapped.Proteins"]<-"Mapped_Proteins"
df_meta$Mapped_Genes<-unique(df_meta$Mapped_Genes)
orfa<-unname(df_meta$MainProtein)
orfa<-gsub(orfa,pattern="R2T2P_",replacement = "")
orfa<-strsplit(orfa,"_")
tx_st<-as.numeric(sapply(orfa,function(x){x[length(x)-1]}))+(df_meta$Start_on_MainProtein-1)*3
tx_orf<-sapply(orfa,"[[",1)
gr_tx_frag<-GRanges(seqnames = tx_orf,ranges = IRanges(start=tx_st,width = nchar(rownames(df_meta))*3))
peps_gen<-pmapFromTranscripts(x = gr_tx_frag,transcripts = GTF_annotation$exons_txs[as.character(seqnames(gr_tx_frag))],ignore.strand=F)
names(peps_gen)<-paste(rownames(df_meta),as.character(gr_tx_frag),sep=";")
peps_gen<-unlist(peps_gen)
peps_gen<-peps_gen[peps_gen$hit]
peps_gen$exon_id<-NULL
peps_gen$exon_rank<-NULL
peps_gen$exon_name<-NULL
peps_gen<-split(peps_gen,names(peps_gen))
names(peps_gen)<-sapply(strsplit(names(peps_gen),";"),"[[",1)
peps_gen<-peps_gen[rownames(df_meta)]
df_meta$coords<-peps_gen[rownames(df_meta)]

# -- Transcript compatibility
ovv<-findOverlaps(df_meta$coords,GTF_annotation$exons_txs,type = "within")
txss<-names(GTF_annotation$exons_txs[ovv@to])
pepss<-rownames(df_meta)[ovv@from]
chunks<-seq(1,length(ovv@from),by = 100000)
if(chunks[length(chunks)]<length(ovv@from)){chunks<-c(chunks,length(ovv@from))}
mapp<-GRangesList()
for(i in 1:(length(chunks)-1)){
  if(i!=(length(chunks)-1)){
    mapp<-suppressWarnings(c(mapp,pmapToTranscripts(df_meta$coords[ovv@from][chunks[i]:(chunks[i+1]-1)],transcripts = GTF_annotation$exons_txs[ovv@to][chunks[i]:(chunks[i+1]-1)])))
  }
  if(i==(length(chunks)-1)){
    mapp<-suppressWarnings(c(mapp,pmapToTranscripts(df_meta$coords[ovv@from][chunks[i]:(chunks[i+1])],transcripts = GTF_annotation$exons_txs[ovv@to][chunks[i]:(chunks[i+1])])))
  }
}
redmapp<-mapp
redmapp<-redmapp[elementNROWS(redmapp)==1]
df_meta$compatible_tx<-CharacterList(unname(split(as.character(seqnames(unlist(redmapp))),names(redmapp))[rownames(df_meta)]))
df_meta$mapping_tx_novelty<-"annotated"
df_meta$mapping_tx_novelty[sum(grepl(df_meta$compatible_tx,pattern="^R1.|^R2."))==elementNROWS(df_meta$compatible_tx)]<-"R1R2"
df_meta$mapping_tx_novelty[sum(grepl(df_meta$compatible_tx,pattern="^R1."))==elementNROWS(df_meta$compatible_tx)]<-"R1"
df_meta$mapping_tx_novelty[sum(grepl(df_meta$compatible_tx,pattern="^R2."))==elementNROWS(df_meta$compatible_tx)]<-"R2"

# -- CDS compatibility
ovv<-findOverlaps(df_meta$coords,GTF_annotation$cds_txs,type = "within")
txss<-names(GTF_annotation$cds_txs[ovv@to])
pepss<-rownames(df_meta)[ovv@from]
chunks<-seq(1,length(ovv@from),by = 100000)
if(chunks[length(chunks)]<length(ovv@from)){chunks<-c(chunks,length(ovv@from))}
mapp<-GRangesList()
for(i in 1:(length(chunks)-1)){
  if(i!=(length(chunks)-1)){
    mapp<-suppressWarnings(c(mapp,pmapToTranscripts(df_meta$coords[ovv@from][chunks[i]:(chunks[i+1]-1)],transcripts = GTF_annotation$cds_txs[ovv@to][chunks[i]:(chunks[i+1]-1)])))
  }
  if(i==(length(chunks)-1)){
    mapp<-suppressWarnings(c(mapp,pmapToTranscripts(df_meta$coords[ovv@from][chunks[i]:(chunks[i+1])],transcripts = GTF_annotation$cds_txs[ovv@to][chunks[i]:(chunks[i+1])])))
  }
}
redmapp<-mapp
redmapp<-redmapp[elementNROWS(redmapp)==1]
compatibility_cds<-CharacterList(split(as.character(seqnames(unlist(redmapp))),names(redmapp)))
nocdspeps<-rownames(df_meta)[!rownames(df_meta)%in%names(compatibility_cds)]
if(length(nocdspeps)>0){
  nocdss<-CharacterList(rep(list(NA),length(nocdspeps)))
  names(nocdss)<-nocdspeps
  compatibility_cds[nocdspeps]<-nocdss
}
df_meta$compatible_cds<-unname(compatibility_cds[rownames(df_meta)])
df_meta$mapping_cds_novelty<-"annotated"
df_meta$mapping_cds_novelty[sum(is.na(df_meta$compatible_cds))>0]<-"novel"

# -- Other
df_meta$Proteins_collapsed<-unlist(lapply(df_meta$Mapped_Proteins,function(x){paste(sort(x),collapse=";")}))
df_meta$gene_id<-unlist(lapply(df_meta$Mapped_Genes,function(x){paste(sort(unique(x)),collapse=";")}))
# Add the PG fragpipe, and go for limma for all and diffsplice peptide/gene, ORFgroup/gene and ProteinGroup/gene
genes_pep<-unname(unique(df_meta$Mapped_Genes))
genes_pep<-unlist(lapply(genes_pep,function(x){paste(x,sep=";")}))
###
# TODO 3 things:
# 1) Add the compatible_biotype column to check if it is coming from R1 or R2 (or both) ORFs
# 2) If only ORFquant entries add NA to the novelty column, but retain R1 or R2 category based on point 1
# 3) Change prioritization in ORFquant code, I know where we have to edit
###
# TODO Change to R2T2P_custom
df_meta$in_R2T2P<-sum(grepl(df_meta$Mapped_Proteins,pattern="R2T2P_"))>0
# TODO Change to R2T2P_ref
df_meta$in_ref<-sum(!grepl(df_meta$Mapped_Proteins,pattern="R2T2P_"))>0
df_meta$novel<-df_meta$in_R2T2P & !df_meta$in_ref
df_meta$in_R2T2P_noveltx<-"none"
df_meta$in_R2T2P_noveltx[sum(grepl(df_meta$compatible_biotype,pattern="novel_tx_R2|novel_tx_R1"))==elementNROWS(df_meta$compatible_biotype) & df_meta$mapping_tx_novelty!="annotated"]<-"novel_tx_R1R2"
df_meta$in_R2T2P_noveltx[sum(grepl(df_meta$compatible_biotype,pattern="novel_tx_R1"))==elementNROWS(df_meta$compatible_biotype) & df_meta$mapping_tx_novelty=="R1"]<-"novel_tx_R1"
df_meta$in_R2T2P_noveltx[sum(grepl(df_meta$compatible_biotype,pattern="novel_tx_R2"))==elementNROWS(df_meta$compatible_biotype) & df_meta$mapping_tx_novelty=="R2"]<-"novel_tx_R2"
df_meta$novel_annotated<-df_meta$novel & df_meta$in_R2T2P_noveltx=="none"
df_meta$novelty_category<-"in_ref"
df_meta$novelty_category[df_meta$novel]<-"novel"
df_meta$novelty_category[df_meta$novel & df_meta$in_R2T2P_noveltx=="novel_tx_R1R2"]<-"novel_tx_R1R2"
df_meta$novelty_category[df_meta$novel & df_meta$in_R2T2P_noveltx=="novel_tx_R1"]<-"novel_tx_R1"
df_meta$novelty_category[df_meta$novel & df_meta$in_R2T2P_noveltx=="novel_tx_R2"]<-"novel_tx_R2"
# Careful here, derive novelty from CDS as well
if(sum(df_meta$in_ref)==0){
  df_meta$novelty_category[df_meta$novel_annotated]<-"annotated_tx"
  df_meta$novelty_category[df_meta$mapping_cds_novelty=="novel" & !grepl(df_meta$novelty_category,pattern="^novel_tx")]<-"novel_CDS"
}
trann<-unique(GTF_annotation$trann[,c("gene_id","gene_name","gene_biotype")])
df_meta$gene_biotype<-trann$gene_biotype[match(df_meta$gene_id,trann$gene_id)]
df_meta$gene_name<-trann$gene_name[match(df_meta$gene_id,trann$gene_id)]
df_meta$gene_biotype[grep(df_meta$gene_biotype,pattern="pseud")]<-"pseudogene"
df_meta$gene_biotype[grep(df_meta$gene_biotype,pattern="novel")]<-"novel"
df_meta$gene_biotype[is.na(df_meta$gene_name)]<-"overlapping"
df_meta$gene_biotype[!df_meta$gene_biotype%in%c("protein_coding","pseudogene","lncRNA","novel","overlapping")]<-"other"
df_meta$gene_id_novelty<-paste(df_meta$gene_id,df_meta$novelty_category,sep=";")

# -- Peptide data (raw intensities + imputed log2-transformed intensities + metadata + imputation fit)
df_int<-peptide_intensities
df_int[df_int == 0] <- NA
df_int_nolog<-df_int
df_int<-apply(df_int,2,log2)
dpcfit <- dpc(df_int)
df_int_imp <- dpcImpute(df_int,dpc = dpcfit,verbose = F)
peptide_data<-list(df_int_nolog,df_int_imp,df_meta,dpcfit)
names(peptide_data)<-c("intensities_raw","intensities_log2_imputed","metadata","imputation_fit")

# -- Peptide GTF file
peps_gtf<-unlist(peptide_data$metadata$coords)
mcols(peps_gtf)<-NULL
pepseq<-sapply(strsplit(names(peps_gtf),"[.]"),"[[",1)
orfid<-sapply(strsplit(names(peps_gtf),";"),"[[",2)
peps_gtf$transcript_id<-pepseq
peps_gtf$gene_id<-orfid
peps_gtf$gene_name<-peptide_data$metadata$gene_name[match(peps_gtf$transcript_id,rownames(peptide_data$metadata))]
peps_gtf$novelty<-peptide_data$metadata$novelty_category[match(peps_gtf$transcript_id,rownames(peptide_data$metadata))]
peps_gtf$type<-"CDS"
peps_gtf$`source`="R2T2P"
names(peps_gtf)<-NULL
# Saving the peptide GTF file
suppressWarnings(export.gff2(object=peps_gtf,con=paste0(search_fold, "_peptide_coords.gtf")))

# -- Additional steps with metadata of peptide data
df_meta<-peptide_data$metadata
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
gns_peps<-split(unstrsplit(unname(df_meta$Mapped_Proteins), sep = ";"),df_meta$gene_id)
multi<-elementNROWS(unique(gns_peps))>1
multi<-gns_peps[multi]
multipeps<-CharacterList(lapply(multi,function(x){
  xx<-CharacterList(strsplit(x,";"))
  unqo<-rep(NA,length(x))
  entries<-NA
  for(i in 1:length(xx)){
    xxx<-xx[[i]]
    unqoxxx<-which(sum(xx%in%xxx)==0)
    if(length(unqoxxx)>0){unqo[i]<-paste(unqoxxx,collapse=";")}
  }
  if(sum(!is.na(unqo))>0){
    unqok<-unqo[!is.na(unqo)]
    unqunqo<-names(sort(table(unqok),decreasing = T))
    unqunqospl<-strsplit(unqunqo,";")
    entries<-unstrsplit(lapply(unqunqospl,function(y){
      unique(unlist(strsplit(x[as.numeric(y)],";")))
    }),";")}
  entries
}))
multiprot<-which(sum(!is.na(multipeps))>0)
multiprot<-multipeps[multiprot]
names(multiprot)<-GTF_annotation$trann$gene_name[match(names(multiprot),GTF_annotation$trann$gene_id)]
top_cand<-sapply(multiprot,"[[",1)
df_meta<-peptide_data$metadata
df_meta$multiprot_candidate<-FALSE
df_meta$multiprot_candidate[unstrsplit(df_meta$Mapped_Proteins,";")%in%top_cand]<-TRUE
peptide_data$metadata<-df_meta
# Saving the peptide data as an RData file
save(peptide_data, file=paste0(search_fold, "_peptide_results_annotated.RData"))

# -- Differential analyses

df_meta<-peptide_data$metadata
multi<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))>1,]
df_int<-peptide_data$intensities_raw
df_int<-log2(as.matrix(df_int))
df_int<-df_int[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]

if (data_type == "DDA_TMT") {
 cond <- de_table$condition
} else {
 cond <- simpl_de_table$condition
}
comparisons<-paste(unique(cond)[unique(cond)!=baseline],"-",baseline,sep="")

#4 Groups: peptides, PG fragpipe, Genes, proteins_collapsed,gene_id_novelty
# for diffsplice peptides|gene, Proteins_collapsed|genes, PGfragpipe|gene
#missing PGfragpipe at the moment, let's go with what we have

prot_groups<-df_meta[,c("gene_id","Proteins_collapsed","gene_id_novelty")]
prot_groups$Peptides<-rownames(df_meta)


DE_res<-list()

for(prot_group in colnames(prot_groups)){
  #make a for look with DE results per each grouping
  if(prot_group=="Peptides"){
    y.protein=as.matrix(df_int)
  }else{
    protein_group<-df_meta[,prot_group]
    #this with imputation
    #y.protein_dpc <- dpcQuant(df_int, protein_group, dpc=dpcfit,verbose = F)
    df_int_new<-2^df_int
    y.protein <- aggregate(df_int_new,by=list(protein_group),sum,na.rm=T)
    rownames(y.protein)<-y.protein$Group.1
    y.protein<-log2(y.protein[,-1])

  }
  design <- model.matrix(~0+cond)
  colnames(design) <- gsub(colnames(design), pattern="cond", replacement="")
  contrast <- makeContrasts(contrasts=comparisons,levels=design)
  #this also with imputation
  #fit <- dpcDE(y.protein, design,plot = F)
  fit = lmFit(y.protein, design)
  fit2 = contrasts.fit(fit, contrast)
  fit2 = eBayes(fit2)

  DE_results <- list()
  for (j in 1:length(comparisons)) {
    #ADD AVG EXPRESSION CALCULATED BY HAND

    DE_result <- topTable(fit2, coef=comparisons[j], number=Inf,sort.by="none")
    DE_result$padj<-DE_result$adj.P.Val
    DE_result<-DE_result[,c("AveExpr","logFC","padj")]
    if(prot_group=="Peptides"){
      coor<-unlist(df_meta$coords)
      coor$logFC<-rep(DE_result[names(df_meta$coords),"logFC"],elementNROWS(df_meta$coords))
      coor$logFC[is.na(coor$logFC)]<-0
      covopl<-coverage(coor[as.character(strand(coor))=="+"],weight = coor[as.character(strand(coor))=="+"]$logFC)
      covomn<-coverage(coor[as.character(strand(coor))=="-"],weight = coor[as.character(strand(coor))=="-"]$logFC)
      runValue(covopl)<-round(runValue(covopl),digits = 4)
      runValue(covomn)<-round(runValue(covomn),digits = 4)
      # Exporting bigwig files
      export.bw(covopl,con = paste0(search_fold,"_peptide_log2FC_",comparisons[j],"_plus.bw"))
      export.bw(covomn,con = paste0(search_fold,"_peptide_log2FC_",comparisons[j],"_minus.bw"))
    }

    DE_results[[j]] <- DE_result
  }
  names(DE_results) <- comparisons
  #if(prot_group!="Peptides"){
  DE_results$intensities<-y.protein
  #}
  DE_res[[prot_group]]<-DE_results
}

#here is about peps, be careful

prot_coll<-unique(df_meta[,c("gene_id","Proteins_collapsed")])
rownames(prot_coll)<-NULL
prot_coll<-prot_coll[order(prot_coll$Proteins_collapsed,decreasing = F),]

#missing PG fragpipe
diffsplice_groups<-list(df_meta$gene_id,prot_coll$gene_id)
names(diffsplice_groups)<-c("Peptides","Proteins_collapsed")
#missing PG fragpipe
df_intensities<-list(DE_res$Peptides$intensities,DE_res$Proteins_collapsed$intensities[unname(prot_coll$Proteins_collapsed),])
names(df_intensities)<-c("Peptides","Proteins_collapsed")
diffsplice_res<-list()
for(gene_group in names(df_intensities)){
  df_vals<-df_intensities[[gene_group]]
  gnns<-unname(diffsplice_groups[[gene_group]])

  design <- model.matrix(~0+cond)
  colnames(design) <- gsub(colnames(design), pattern="cond", replacement="")
  contrast <- makeContrasts(contrasts=comparisons,levels=design)

  fit = lmFit(df_vals, design)
  fit = contrasts.fit(fit, contrast)
  fit = eBayes(fit)

  diff_results <- list()
  for (jj in 1:length(comparisons)) {
    # Getting DE results + adding DE status

    res_diffsplice <- diffSplice(fit, geneid=gnns,verbose = F,robust = F)
    res_diffsplice <- topSplice(res_diffsplice, coef=comparisons[jj], number=Inf,test = "t",sort.by="none")
    res_diffsplice$padj<-res_diffsplice$FDR
    #this needs checking
    #if(gene_group=="Peptides"){
    rownames(res_diffsplice)<-rownames(df_vals)[as.numeric(rownames(res_diffsplice))]
    #}
    df_g<-split(DE_res[[gene_group]][[jj]],f=gnns)
    df_g<-lapply(df_g,function(x){
      fcs<-x$logFC
      deltafcs<-c()
      for(j in 1:length(fcs)){
        meanfcs<-mean(fcs[-j],na.rm=T)
        deltafcs<-c(deltafcs,fcs[j]-meanfcs)
      }
      x$deltafcs<-deltafcs
      x$meanfcs<-meanfcs
      x$id<-rownames(x)
      rownames(x)<-NULL
      x
    })
    df_g<-do.call(df_g,what=rbind)
    rownames(df_g)<-df_g$id
    res_diffsplice$delta_log2FC<-df_g[rownames(res_diffsplice),]$deltafcs
    res_diffsplice$other_log2FC<-df_g[rownames(res_diffsplice),]$meanfcs

    ints_g<-DE_res[[gene_group]]$intensities
    cond1<-strsplit(comparisons[jj],"-")[[1]][1]
    cond2<-strsplit(comparisons[jj],"-")[[1]][2]

    df_g<-split(data.frame(ints_g),f=gnns)
    df_g<-lapply(df_g,function(x){
      xx<-as.matrix(x)
      fcs<-rowMeans(xx[,design[,cond1]>0,drop=F],na.rm=T)-rowMeans(xx[,design[,cond2]>0,drop=F],na.rm=T)
      deltafcs<-c()
      for(j in 1:length(fcs)){
        meanfcs<-mean(fcs[-j],na.rm=T)
        deltafcs<-c(deltafcs,fcs[j]-meanfcs)
      }
      x$deltafcs<-deltafcs
      x$meanfcs<-meanfcs
      x$id<-rownames(x)
      rownames(x)<-NULL
      x
    })
    df_g<-do.call(df_g,what=rbind)
    rownames(df_g)<-df_g$id
    res_diffsplice$delta_log2FC_raw<-df_g[rownames(res_diffsplice),]$deltafcs
    res_diffsplice$other_log2FC_raw<-df_g[rownames(res_diffsplice),]$meanfcs

    res_diffsplice$group_id<-res_diffsplice$GeneID
    res_diffsplice<-res_diffsplice[,c("group_id","logFC","padj","delta_log2FC","other_log2FC")]

    missingpeps<-rownames(df_vals)[!rownames(df_vals)%in%rownames(res_diffsplice)]
    diffres<-matrix(NA,nrow = length(missingpeps),ncol = ncol(res_diffsplice))
    diffres<-data.frame(diffres)
    colnames(diffres)<-colnames(res_diffsplice)
    rownames(diffres)<-missingpeps
    res_diffsplice<-rbind(res_diffsplice,diffres)[rownames(df_vals),]

    diff_results[[jj]] <- res_diffsplice
  }
  names(diff_results) <- comparisons
  diffsplice_res[[gene_group]]<-diff_results
}

#add the not-imputed results

list_DE_raw<-list(DE_res,diffsplice_res)
names(list_DE_raw)<-c("DE","diffsplice")

df_meta<-peptide_data$metadata
multi<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))>1,]
df_int<-peptide_data$intensities_log2_imputed$E
#df_int<-log2(as.matrix(df_int))

df_int<-df_int[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]


#make tables, then go for DE, then group also by DUX4 target category (from gene ids)
#supplementary for imputed values


#now group them into protein groups, according to protein groups (fragpipe-derived),
#genes, genes dividing between novel peptides and novel peptides in novel txs
#and isoform groups (not sure how though...maybe get inspired by ORFquant)

if (data_type == "DDA_TMT") {
  cond <- de_table$condition
} else {
  cond <- simpl_de_table$condition
}

#4 Groups: peptides, PG fragpipe, Genes, proteins_collapsed,gene_id_novelty
# for diffsplice peptides|gene, Proteins_collapsed|genes, PGfragpipe|gene

#missing PGfragpipe at the moment, let's go with what we have

prot_groups<-df_meta[,c("gene_id","Proteins_collapsed","gene_id_novelty")]
prot_groups$Peptides<-rownames(df_meta)

DE_res<-list()

for(prot_group in colnames(prot_groups)){
  #make a for look with DE results per each grouping
  if(prot_group=="Peptides"){
    y.protein=as.matrix(df_int)
  }else{
    protein_group<-df_meta[,prot_group]
    #this with imputation
    #y.protein_dpc <- dpcQuant(df_int, protein_group, dpc=dpcfit,verbose = F)
    df_int_new<-2^df_int
    y.protein <- aggregate(df_int_new,by=list(protein_group),sum,na.rm=T)
    rownames(y.protein)<-y.protein$Group.1
    y.protein<-log2(y.protein[,-1])

  }
  design <- model.matrix(~0+cond)
  colnames(design) <- gsub(colnames(design), pattern="cond", replacement="")
  contrast <- makeContrasts(contrasts=comparisons,levels=design)
  #this also with imputation
  #fit <- dpcDE(y.protein, design,plot = F)
  fit = lmFit(y.protein, design)
  fit2 = contrasts.fit(fit, contrast)
  fit2 = eBayes(fit2)

  DE_results <- list()
  for (j in 1:length(comparisons)) {
    # Getting DE results + adding DE status
    DE_result <- topTable(fit2, coef=comparisons[j], number=Inf,sort.by="none")
    DE_result$padj<-DE_result$adj.P.Val
    DE_result<-DE_result[,c("AveExpr","logFC","padj")]
    DE_results[[j]] <- DE_result
  }
  names(DE_results) <- comparisons
  DE_results$intensities<-y.protein
  DE_res[[prot_group]]<-DE_results
}

#here is about peps, be careful

prot_coll<-unique(df_meta[,c("gene_id","Proteins_collapsed")])
rownames(prot_coll)<-NULL
prot_coll<-prot_coll[order(prot_coll$Proteins_collapsed,decreasing = F),]

#missing PG fragpipe
diffsplice_groups<-list(df_meta$gene_id,prot_coll$gene_id)
names(diffsplice_groups)<-c("Peptides","Proteins_collapsed")
#missing PG fragpipe
df_intensities<-list(DE_res$Peptides$intensities,DE_res$Proteins_collapsed$intensities[unname(prot_coll$Proteins_collapsed),])
names(df_intensities)<-c("Peptides","Proteins_collapsed")
diffsplice_res<-list()
for(gene_group in names(df_intensities)){
  df_vals<-df_intensities[[gene_group]]
  gnns<-unname(diffsplice_groups[[gene_group]])


  design <- model.matrix(~0+cond)
  colnames(design) <- gsub(colnames(design), pattern="cond", replacement="")
  contrast <- makeContrasts(contrasts=comparisons,levels=design)

  fit = lmFit(df_vals, design)
  fit = contrasts.fit(fit, contrast)
  fit = eBayes(fit)

  diff_results <- list()
  for (jj in 1:length(comparisons)) {
    # Getting DE results + adding DE status

    res_diffsplice <- diffSplice(fit, geneid=gnns,verbose = F,robust = F)
    res_diffsplice <- topSplice(res_diffsplice, coef=comparisons[jj], number=Inf,test = "t",sort.by="none")
    res_diffsplice$padj<-res_diffsplice$FDR
    #this needs checking
    #if(gene_group=="Peptides"){
    rownames(res_diffsplice)<-rownames(df_vals)[as.numeric(rownames(res_diffsplice))]
    #}
    df_g<-split(DE_res[[gene_group]][[jj]],f=gnns)
    df_g<-lapply(df_g,function(x){
      fcs<-x$logFC
      deltafcs<-c()
      for(j in 1:length(fcs)){
        meanfcs<-mean(fcs[-j],na.rm=T)
        deltafcs<-c(deltafcs,fcs[j]-meanfcs)
      }
      x$deltafcs<-deltafcs
      x$meanfcs<-meanfcs
      x$id<-rownames(x)
      rownames(x)<-NULL
      x
    })
    df_g<-do.call(df_g,what=rbind)
    rownames(df_g)<-df_g$id
    res_diffsplice$delta_log2FC<-df_g[rownames(res_diffsplice),]$deltafcs
    res_diffsplice$other_log2FC<-df_g[rownames(res_diffsplice),]$meanfcs

    ints_g<-DE_res[[gene_group]]$intensities
    cond1<-strsplit(comparisons[jj],"-")[[1]][1]
    cond2<-strsplit(comparisons[jj],"-")[[1]][2]

    df_g<-split(data.frame(ints_g),f=gnns)
    df_g<-lapply(df_g,function(x){
      xx<-as.matrix(x)
      fcs<-rowMeans(xx[,design[,cond1]>0,drop=F],na.rm=T)-rowMeans(xx[,design[,cond2]>0,drop=F],na.rm=T)
      deltafcs<-c()
      for(j in 1:length(fcs)){
        meanfcs<-mean(fcs[-j],na.rm=T)
        deltafcs<-c(deltafcs,fcs[j]-meanfcs)
      }
      x$deltafcs<-deltafcs
      x$meanfcs<-meanfcs
      x$id<-rownames(x)
      rownames(x)<-NULL
      x
    })
    df_g<-do.call(df_g,what=rbind)
    rownames(df_g)<-df_g$id
    res_diffsplice$delta_log2FC_raw<-df_g[rownames(res_diffsplice),]$deltafcs
    res_diffsplice$other_log2FC_raw<-df_g[rownames(res_diffsplice),]$meanfcs


    res_diffsplice$group_id<-res_diffsplice$GeneID
    res_diffsplice<-res_diffsplice[,c("group_id","logFC","padj","delta_log2FC","other_log2FC","delta_log2FC_raw","other_log2FC_raw")]

    missingpeps<-rownames(df_vals)[!rownames(df_vals)%in%rownames(res_diffsplice)]
    diffres<-matrix(NA,nrow = length(missingpeps),ncol = ncol(res_diffsplice))
    diffres<-data.frame(diffres)
    colnames(diffres)<-colnames(res_diffsplice)
    rownames(diffres)<-missingpeps
    res_diffsplice<-rbind(res_diffsplice,diffres)[rownames(df_vals),]

    diff_results[[jj]] <- res_diffsplice
  }
  names(diff_results) <- comparisons
  diffsplice_res[[gene_group]]<-diff_results
}

list_DE_imp<-list(DE_res,diffsplice_res)
names(list_DE_imp)<-c("DE","diffsplice")

df_meta<-peptide_data$metadata
multi<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))>1,]
df_int<-peptide_data$intensities_log2_imputed
#df_int<-log2(as.matrix(df_int))

df_int<-df_int[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]


#make tables, then go for DE, then group also by DUX4 target category (from gene ids)
#supplementary for imputed values


#now group them into protein groups, according to protein groups (fragpipe-derived),
#genes, genes dividing between novel peptides and novel peptides in novel txs
#and isoform groups (not sure how though...maybe get inspired by ORFquant)

if (data_type == "DDA_TMT") {
  cond <- de_table$condition
} else {
  cond <- simpl_de_table$condition
}

#4 Groups: peptides, PG fragpipe, Genes, proteins_collapsed,gene_id_novelty
# for diffsplice peptides|gene, Proteins_collapsed|genes, PGfragpipe|gene

#missing PGfragpipe at the moment, let's go with what we have

prot_groups<-df_meta[,c("gene_id","Proteins_collapsed","gene_id_novelty")]
prot_groups$Peptides<-rownames(df_meta)

DE_res<-list()

for(prot_group in colnames(prot_groups)){
  #make a for look with DE results per each grouping
  if(prot_group=="Peptides"){
    y.protein=df_int
  }else{
    protein_group<-df_meta[,prot_group]
    #this with imputation
    y.protein <- dpcQuant(df_int, protein_group, dpc=peptide_data$imputation_fit,verbose = F)
    #df_int_new<-2^df_int
    #y.protein <- aggregate(df_int_new,by=list(protein_group),sum,na.rm=T)
    #rownames(y.protein)<-y.protein$Group.1
    #y.protein<-log2(y.protein[,-1])

  }
  design <- model.matrix(~0+cond)
  colnames(design) <- gsub(colnames(design), pattern="cond", replacement="")
  contrast <- makeContrasts(contrasts=comparisons,levels=design)
  #this also with imputation
  #fit <- dpcDE(y.protein, design,plot = F)
  fit = dpcDE(y.protein, design, plot=F)
  fit2 = contrasts.fit(fit, contrast)
  fit2 = eBayes(fit2)

  DE_results <- list()
  for (j in 1:length(comparisons)) {
    # Getting DE results + adding DE status
    DE_result <- topTable(fit2, coef=comparisons[j], number=Inf,sort.by="none")
    DE_result$padj<-DE_result$adj.P.Val
    DE_result<-DE_result[,c("AveExpr","logFC","padj")]
    DE_results[[j]] <- DE_result
  }
  names(DE_results) <- comparisons
  #if(prot_group!="Peptides"){
  DE_results$intensities<-y.protein
  #}
  DE_res[[prot_group]]<-DE_results
}

#here is about peps, be careful

prot_coll<-unique(df_meta[,c("gene_id","Proteins_collapsed")])
rownames(prot_coll)<-NULL
prot_coll<-prot_coll[order(prot_coll$Proteins_collapsed,decreasing = F),]

#missing PG fragpipe
diffsplice_groups<-list(df_meta$gene_id,prot_coll$gene_id)
names(diffsplice_groups)<-c("Peptides","Proteins_collapsed")
#missing PG fragpipe
df_intensities<-list(DE_res$Peptides$intensities,DE_res$Proteins_collapsed$intensities[unname(prot_coll$Proteins_collapsed),])
names(df_intensities)<-c("Peptides","Proteins_collapsed")
diffsplice_res<-list()
for(gene_group in names(df_intensities)){
  df_vals<-df_intensities[[gene_group]]
  gnns<-unname(diffsplice_groups[[gene_group]])

  design <- model.matrix(~0+cond)
  colnames(design) <- gsub(colnames(design), pattern="cond", replacement="")
  contrast <- makeContrasts(contrasts=comparisons,levels=design)
  if(gene_group=="Peptides"){
    fit = dpcDE(df_int, design, plot=F)
    df_vals<-df_vals$E
  }else{
    fit = lmFit(df_vals, design)
  }
  fit = contrasts.fit(fit, contrast)
  fit = eBayes(fit)

  diff_results <- list()
  for (jj in 1:length(comparisons)) {
    # Getting DE results + adding DE status

    res_diffsplice <- diffSplice(fit, geneid=gnns,verbose = F,robust = F)
    res_diffsplice <- topSplice(res_diffsplice, coef=comparisons[jj], number=Inf,test = "t",sort.by="none")
    res_diffsplice$padj<-res_diffsplice$FDR
    #this needs checking
    if(gene_group=="Peptides"){
      rownames(res_diffsplice)<-rownames(peptide_data$intensities_log2_imputed$E)[as.numeric(rownames(res_diffsplice))]
      #add NA entries for peps with no diffsplice (one pep per gene)
      xxx<-res_diffsplice[1,]
      for(jjj in 1:dim(xxx)[2]){
        xxx[,jjj]<-NA
      }
      rnms_miss<-rownames(df_vals)[!rownames(df_vals)%in%rownames(res_diffsplice)]
      xxx[1:length(rnms_miss),]<-xxx
      rownames(xxx)<-rnms_miss
      res_diffsplice<-rbind(res_diffsplice,xxx)

    }

    df_g<-split(DE_res[[gene_group]][[jj]],f=gnns)
    df_g<-lapply(df_g,function(x){
      fcs<-x$logFC
      deltafcs<-c()
      for(j in 1:length(fcs)){
        meanfcs<-mean(fcs[-j],na.rm=T)
        deltafcs<-c(deltafcs,fcs[j]-meanfcs)
      }
      x$deltafcs<-deltafcs
      x$meanfcs<-meanfcs
      x$id<-rownames(x)
      rownames(x)<-NULL
      x
    })
    df_g<-do.call(df_g,what=rbind)
    rownames(df_g)<-df_g$id
    res_diffsplice$delta_log2FC<-df_g[rownames(res_diffsplice),]$deltafcs
    res_diffsplice$other_log2FC<-df_g[rownames(res_diffsplice),]$meanfcs

    ints_g<-DE_res[[gene_group]]$intensities
    cond1<-strsplit(comparisons[jj],"-")[[1]][1]
    cond2<-strsplit(comparisons[jj],"-")[[1]][2]

    df_g<-split(data.frame(ints_g),f=gnns)
    df_g<-lapply(df_g,function(x){
      xx<-as.matrix(x)
      fcs<-rowMeans(xx[,design[,cond1]>0,drop=F],na.rm=T)-rowMeans(xx[,design[,cond2]>0,drop=F],na.rm=T)
      deltafcs<-c()
      for(j in 1:length(fcs)){
        meanfcs<-mean(fcs[-j],na.rm=T)
        deltafcs<-c(deltafcs,fcs[j]-meanfcs)
      }
      x$deltafcs<-deltafcs
      x$meanfcs<-meanfcs
      x$id<-rownames(x)
      rownames(x)<-NULL
      x
    })
    df_g<-do.call(df_g,what=rbind)
    rownames(df_g)<-df_g$id
    res_diffsplice$delta_log2FC_raw<-df_g[rownames(res_diffsplice),]$deltafcs
    res_diffsplice$other_log2FC_raw<-df_g[rownames(res_diffsplice),]$meanfcs
    #took this away to avoid NAs
    #res_diffsplice<-res_diffsplice[rownames(df_vals),]
    res_diffsplice$group_id<-res_diffsplice$GeneID
    res_diffsplice<-res_diffsplice[,c("group_id","logFC","padj","delta_log2FC","other_log2FC","delta_log2FC_raw","other_log2FC_raw")]

    missingpeps<-rownames(df_vals)[!rownames(df_vals)%in%rownames(res_diffsplice)]
    diffres<-matrix(NA,nrow = length(missingpeps),ncol = ncol(res_diffsplice))
    diffres<-data.frame(diffres)
    colnames(diffres)<-colnames(res_diffsplice)
    rownames(diffres)<-missingpeps
    res_diffsplice<-rbind(res_diffsplice,diffres)[rownames(df_vals),]

    diff_results[[jj]] <- res_diffsplice
  }
  names(diff_results) <- comparisons
  diffsplice_res[[gene_group]]<-diff_results
}

list_DE_dcp<-list(DE_res,diffsplice_res)
names(list_DE_dcp)<-c("DE","diffsplice")
list_DE<-list(list_DE_raw,list_DE_imp,list_DE_dcp)
names(list_DE)<-c("raw","imputed","dcp_limpa")

list_intensities<-lapply(list_DE,function(x){
  intens<-lapply(x$DE,"[[","intensities")
  intens$Peptides<-data.frame(intens$Peptides)
  intens
})

list_DE<-lapply(list_DE,function(x){
  x$DE<-lapply(x$DE,function(y){
    y$intensities<-NULL
    y
  })
  x
})

# Saving list_intensities and list_DE as RData files
save(list_intensities,file = paste0(search_fold, "_proteomics_all_intensities.RData"))
save(list_DE,file = paste0(search_fold, "_proteomics_multiDE.RData"))
