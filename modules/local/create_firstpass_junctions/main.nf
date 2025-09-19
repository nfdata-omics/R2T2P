process CREATE_FIRSTPASS_JUNCTIONS {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/nfdata/riboseqc:v1.0.0-patched"

    input:
    tuple val(meta), path(pass1_splice_junctions)
    path "R-user-lib/*"
    path gtf_Rannotation

    output:
    tuple val(meta), path("firstpass_junctions.txt"), emit: pass1_junctions
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export R_LIBS_USER=\$PWD/R-user-lib

    Rscript --vanilla -e '
        library(RiboseQC)

        # Loading RiboseQC annotation
        load_annotation("${gtf_Rannotation}")

        # Getting paths of "SJ.out.tab" files
        sj_out_files_firstpass <- list("${pass1_splice_junctions}")

        # Reading files and storing the junctions (as GRanges objects) inside a list
        firstpass_junctions <- lapply(sj_out_files_firstpass, readSTARJunctions)
        # Converting the list to a GRangesList
        firstpass_junctions <- GRangesList(firstpass_junctions)
        # Unlisting the GRangesList
        firstpass_junctions <- unlist(firstpass_junctions)
        # Removing duplicates
        firstpass_junctions <- unique(firstpass_junctions)

        # Filtering to get truly novel junctions
        ann<-GTF_annotation\$junctions
        nov<-firstpass_junctions
        mcols(ann)<-NULL
        mcols(nov)<-NULL
        firstpass_junctions<-firstpass_junctions[!nov%in%ann]

        # Creating a dataframe with firstpass junctions
        firstpass_junctions_df <- data.frame(seqnames=as.character(seqnames(firstpass_junctions)),
                                            start=as.character(start(firstpass_junctions)),
                                            end=as.character(end(firstpass_junctions)),
                                            strand=as.character(strand(firstpass_junctions)))
        # Replacing values of the column strand according to the manual of STAR 2.7.9a
        firstpass_junctions_df\$strand[firstpass_junctions_df\$strand=="*"]=0
        firstpass_junctions_df\$strand[firstpass_junctions_df\$strand=="+"]=1
        firstpass_junctions_df\$strand[firstpass_junctions_df\$strand=="-"]=2

        # Creating a txt file with novel first-pass junctions
        write.table(firstpass_junctions_df, file="firstpass_junctions.txt", quote=FALSE, sep="\\t", row.names=F, col.names=F)

        # Writing package versions to versions.yml
        x = sessionInfo()
        versions <- list(
            R = paste0(x\$R.version\$major, ".", x\$R.version\$minor)
        )
        for(i in seq_along(x\$otherPkgs)){
            pkg <- x\$otherPkgs[[i]]
            versions[[pkg\$Package]] <- pkg\$Version
        }
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

    stub:
    """
    touch firstpass_junctions.txt

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
        # Convert list to yaml and write to file
        yaml::write_yaml(versions, "versions.yml")
    '
    """

}
