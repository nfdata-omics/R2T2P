process CREATE_PROTEIN_DB_WRITE_DB {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/nfdata/riboseqc:v1.3.0-patched"

    input:
    path "R-user-lib/*"
    path gtf_Rannot
    path gtf_stringtie_Rannot
    tuple val(meta), path(orfquant_protein_fasta)

    output:
    path "annot_proteins_db.fasta",              emit: annot_prot_fasta
    path "orfquant_proteins_db.fasta",           emit: orfquant_prot_fasta
    path "annot_and_orfquant_proteins_db.fasta", emit: annot_and_orfquant_prot_fasta
    path "prot_ID_pairs.txt",                    emit: prot_id_pairs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export R_LIBS_USER=\$PWD/R-user-lib

    Rscript --vanilla -e '
        # Script to create the protein database (annotated proteins + ORFquant proteins)

        # --- FASTA file of annotated proteins

        # Attaching the R package RiboseQC
        suppressMessages(library("RiboseQC"))

        # Loading the Rannot
        load_annotation("${gtf_Rannot}")

        # Getting the protein sequences
        seqs_coding <- extractTranscriptSeqs(x = genome_seq, transcripts = GTF_annotation\$cds_txs)
        annot_proteins <- suppressWarnings(translate(seqs_coding, if.fuzzy.codon = "solve"))

        chrs <- as.vector(seqnames(unlist(range(GTF_annotation\$cds_txs))))
        mitos <- which(as.data.frame(seqinfo(GTF_annotation\$cds_txs))[chrs,"isCircular"])

        annot_proteins[mitos] <- suppressWarnings(translate(seqs_coding[mitos], genetic.code = getGeneticCode("2"), if.fuzzy.codon = "solve"))
        cdsss <- unname(GTF_annotation\$cds_txs_coords)
        mcols(cdsss) <- NULL
        nmss <- as.character(unstrand(cdsss))
        nmss <- paste(CharacterList(strsplit(nmss,"[:-]")), collapse ="_")

        # Removing all asterisks at the end of protein sequences
        last_chars <- as.character(subseq(annot_proteins, start=width(annot_proteins), end=width(annot_proteins)))
        while ("*" %in% last_chars) {
            annot_proteins[last_chars == "*"] <- subseq(annot_proteins[last_chars == "*"], start=1, end=width(annot_proteins[last_chars == "*"])-1)
            last_chars <- as.character(subseq(annot_proteins, start=width(annot_proteins), end=width(annot_proteins)))
        }

        # Adding names of protein sequences
        names(annot_proteins) <- nmss

        # --- FASTA file of ORFquant proteins

        orfquant_proteins <- readAAStringSet("${orfquant_protein_fasta}", format="fasta")
        names(orfquant_proteins) <- sapply(strsplit(names(orfquant_proteins), split="\\\\|"), "[[", 1)
        names(orfquant_proteins) <- paste("R2T2P_", names(orfquant_proteins), sep="")

        # --- FASTA file of annotated and ORFquant proteins

        # Combining the two AAStringSet objects (annotated proteins + ORFquant proteins)
        annot_and_orfquant_proteins <- c(annot_proteins, orfquant_proteins)

        # --- Replace original FASTA headers with Uniprot-like FASTA headers

        # Loading the Rannot (annotated + de novo transcripts)
        load_annotation("${gtf_stringtie_Rannot}")

        # Defining new protein IDs and creating a txt file with pairs (original and new protein IDs)
        new_prot_ids <- paste0("P", 100000000 + 1:length(annot_and_orfquant_proteins))
        prot_ID_pairs_df <- data.frame(orig_id=names(annot_and_orfquant_proteins), new_id=new_prot_ids)
        write.table(prot_ID_pairs_df, file="prot_ID_pairs.txt", sep="\t", row.names=F, col.names=T, quote=F)

        # Replacing original FASTA headers with Uniprot-like FASTA headers
        new_names_part1 <- paste("tr", new_prot_ids, paste0(new_prot_ids,"_HUMAN"), sep="|")
        new_names_part2 <- "Protein description OS=Homo sapiens OX=9606"
        old_names <- names(annot_and_orfquant_proteins)
        old_names_no_pattern <- gsub(old_names, pattern="R2T2P_", replacement="")
        new_names_part3 <- paste0("GN=", GTF_annotation\$trann\$gene_id[match(sapply(strsplit(old_names_no_pattern, split="_"), "[[", 1), GTF_annotation\$trann\$transcript_id)])
        new_names_part4 <- "PE=1 SV=1"
        new_names <- paste(new_names_part1,
                        new_names_part2,
                        new_names_part3,
                        new_names_part4)
        new_names[grepl(old_names, pattern="^R2T2P")] <- gsub(new_names[grepl(old_names, pattern="^R2T2P")],
                                                            pattern="PE=1",
                                                            replacement="PE=4")
        names(annot_and_orfquant_proteins) <- new_names

        # Filtering to get proteins of the two smaller databases
        annot_proteins <- annot_and_orfquant_proteins[old_names %in% names(annot_proteins)]
        orfquant_proteins <- annot_and_orfquant_proteins[old_names %in% names(orfquant_proteins)]

        # Creating FASTA files of the 3 databases (annotated + ORFquant proteins, annotated proteins, ORFquant proteins)
        writeXStringSet(annot_and_orfquant_proteins, filepath="annot_and_orfquant_proteins_db.fasta", format="fasta")
        writeXStringSet(annot_proteins, filepath="annot_proteins_db.fasta", format="fasta")
        writeXStringSet(orfquant_proteins, filepath="orfquant_proteins_db.fasta", format="fasta")

        # Writing package versions to versions.yml
        x = sessionInfo()
        versions <- list(
            R = paste0(x\$R.version\$major, ".", x\$R.version\$minor)
        )
        for(i in seq_along(x\$otherPkgs)){
            pkg <- x\$otherPkgs[[i]]
            versions[[pkg\$Package]] <- pkg\$Version
        }
        versions <- list("${task.process}" = versions)
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

    stub:
    """
    touch annot_proteins_db.fasta
    touch orfquant_proteins_db.fasta
    touch annot_and_orfquant_proteins_db.fasta
    touch prot_ID_pairs.txt

    Rscript -e '
        library(RiboseQC)

        # Writing package versions to versions.yml
        x = sessionInfo()
        versions <- list(
            R = paste0(x\$R.version\$major, ".", x\$R.version\$minor)
        )
        for(i in seq_along(x\$otherPkgs)){
            pkg <- x\$otherPkgs[[i]]
            versions[[pkg$Package]] <- pkg\$Version
        }
        versions <- list("${task.process}" = versions)
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

}
