process CREATE_PROTEIN_DB_WRITE_DB {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/nfdata/riboseqc:v1.3.0-patched"

    input:
    path "R-user-lib/*"
    path gtf_Rannot
    tuple val(meta), path(orfquant_protein_fasta)

    output:
    path "annot_proteins_db.fasta",              emit: annot_prot_fasta
    path "orfquant_proteins_db.fasta",           emit: orfquant_prot_fasta
    path "annot_and_orfquant_proteins_db.fasta", emit: annot_and_orfquant_prot_fasta
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

        # Creating the FASTA file of annotated proteins
        writeXStringSet(annot_proteins, filepath="annot_proteins_db.fasta", format="fasta")

        # --- FASTA file of ORFquant proteins

        orfquant_proteins <- readAAStringSet("${orfquant_protein_fasta}", format="fasta")
        names(orfquant_proteins) <- sapply(strsplit(names(orfquant_proteins), split="\\\\|"), "[[", 1)
        names(orfquant_proteins) <- paste("R2T2P_", names(orfquant_proteins), sep="")

        # Creating the FASTA file of ORFquant proteins
        writeXStringSet(orfquant_proteins, filepath="orfquant_proteins_db.fasta", format="fasta")

        # --- FASTA file of annotated and ORFquant proteins

        # Combining the two AAStringSet objects (annotated proteins + ORFquant proteins)
        annot_and_orfquant_proteins <- c(annot_proteins, orfquant_proteins)

        # Creating the FASTA file of annotated and ORFquant proteins
        writeXStringSet(annot_and_orfquant_proteins, filepath="annot_and_orfquant_proteins_db.fasta", format="fasta")

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
