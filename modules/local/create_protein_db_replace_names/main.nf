process CREATE_PROTEIN_DB_REPLACE_NAMES {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/nfdata/riboseqc:v1.3.0-patched"

    input:
    tuple val(meta), path("db_fasta.fasta"), path("fasta_with_decoys_and_contam.fasta")

    output:
    tuple val(meta), path("decoys-contam-*.fasta.fas"), emit: fasta_with_decoys_and_contam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export R_LIBS_USER=\$PWD/R-user-lib

    Rscript --vanilla -e '
        # --- Script to replace names of contaminants and reversed contaminants

        # Attaching the R package Biostrings
        suppressMessages(library("Biostrings"))

        # Storing paths to the starting database and to the database with decoys and contaminants
        db_path <- "db_fasta.fasta"
        db_with_decoys_and_contams_path <- "fasta_with_decoys_and_contam.fasta"

        # Reading the starting database and the database with decoys and contaminants
        db <- readAAStringSet(db_path, format="fasta")
        db_with_decoys_and_contams <- readAAStringSet(db_with_decoys_and_contams_path, format="fasta")

        # Getting sequences of contaminants and reversed contaminants
        decoys <- names(db_with_decoys_and_contams)[grepl(names(db_with_decoys_and_contams), pattern="^rev_")]
        contams <- names(db_with_decoys_and_contams)[!names(db_with_decoys_and_contams) %in% c(decoys, names(db))]
        rev_contams <- decoys[gsub(decoys, pattern="^rev_", replacement="") %in% contams]
        contams_and_rev_contams <- db_with_decoys_and_contams[names(db_with_decoys_and_contams) %in% c(contams, rev_contams)]

        # Removing sequences of contaminants and reversed contaminants from the database
        db_without_contams_and_rev_contams <- db_with_decoys_and_contams[!names(db_with_decoys_and_contams) %in% names(contams_and_rev_contams)]

        # Replacing names of contaminants and reversed contaminants
        starting_with_rev <- grepl(names(contams_and_rev_contams), pattern="^rev_")
        names(contams_and_rev_contams)[starting_with_rev] <- gsub(names(contams_and_rev_contams)[starting_with_rev], pattern="^rev_", replacement="rev_contam_")
        names(contams_and_rev_contams)[!starting_with_rev] <- paste0("contam_", names(contams_and_rev_contams)[!starting_with_rev])

        # Adding contaminants and reversed contaminants (with new names) to the database
        db_with_decoys_and_contams <- c(db_without_contams_and_rev_contams, contams_and_rev_contams)

        # Creating the FASTA file
        writeXStringSet(db_with_decoys_and_contams, filepath="decoys-contam-${meta.id}.fasta.fas", format="fasta")

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
    touch renamed_${fasta_with_decoys_and_contam}

    Rscript -e '
        # Attaching the R package Biostrings
        suppressMessages(library("Biostrings"))

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
