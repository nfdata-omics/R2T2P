process FRAGPIPE {
    tag "${meta.id}"
    label 'process_high', 'process_high_memory'

    container "docker.io/fcyucn/fragpipe:24.0"

    input:
    tuple val(meta), path(protein_db_fasta)
    path workflow_file
    path manifest_file
    path "all_mzml/*"
    path tmt_annotation_file
    path tools_folder
    path diann_folder

    output:
    path "*.pepindex",   emit: pepindex
    path "*_search",    emit: fragpipe_results
    path "*.workflow",   emit: workflow_file
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # prepare fragpipe workflow file with the correct database path
    sed "s|^database\\.db-path=.*\$|database.db-path=${protein_db_fasta}|" ${workflow_file} > ${meta.id}.workflow

    # prepare manifest file with the mzML paths where links will be created organized by experiment_name
    awk -F'\t' 'BEGIN{OFS=FS} {n=split(\$1,a,"/"); d=\$2; \$1=ENVIRON["PWD"] "/" d "/" a[n]; print}' ${manifest_file} > manifest_file_new_paths

    # create experiment_name directories and add to each subfolder the annotation file
    for exp_name in \$(awk '{print \$2}' manifest_file_new_paths  | sort -u); do
        mkdir -p "\$exp_name"
        if [ "\$(basename "${tmt_annotation_file}")" != "NO_FILE" ]; then
            awk -v s="_\$exp_name" '{ if (\$2 != "NA") \$2 = \$2 s; print \$1, \$2 }' "$tmt_annotation_file" > "\$exp_name/tmt_annotation.txt"
        fi
    done

    # create links to the staged mzML files in each experiment subfolder
    awk '{n=split(\$1,a,"/"); fname=a[n]; print "cp -L \\"all_mzml/" fname "\\" \\"" \$1 "\\"" }' manifest_file_new_paths | bash

    # manually stage with copy the external tools and diann installation
    cp -Lr ${tools_folder} ./tools-copy
    cp -Lr ${diann_folder} ./diann-copy

    # lauch fragpipe
    FRAGPIPE_EXE="/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/bin/fragpipe"
    \$FRAGPIPE_EXE \
                --headless \
                --workflow ${meta.id}.workflow \
                --manifest manifest_file_new_paths \
                --workdir ${meta.id}_search \
                --config-tools-folder ./tools-copy \
                --config-diann ./diann-copy/diann-linux \
                --config-python /usr/bin/python3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fragpipe: \$(\${FRAGPIPE_EXE} -h | head -n 1 | awk '{print \$2}' )
        diann: \$(${diann_folder}/diann-linux --version |  --version | awk '(NR==2) {print \$2}' )
        IonQuant: \$(java -jar ${tools_folder}/IonQuant-*/IonQuant-*.jar | head -n1 | awk -F- '{print \$2}')
        MSFragger: \$(java -jar ${tools_folder}/MSFragger-*/MSFragger-*.jar | head -n1 | awk -F- '{print \$2}')
        diaTracer: \$(java -jar ${tools_folder}/diaTracer-*/diaTracer-*.jar | head -n1 | awk -F- '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    mkdir fragpipe_${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fragpipe: \$(\${FRAGPIPE_EXE} -h | head -n 1 | awk '{print \$2}' )
        diann: \$(\${diann_folder}/diann-linux --version |  --version | awk '(NR==2) {print \$2}' )
        IonQuant: \$(java -jar ${tools_folder}/IonQuant-*.jar | head -n1 | awk -F- '{print \$2}')
        MSFragger: \$(java -jar ${tools_folder}/MSFragger-*.jar | head -n1 | awk -F- '{print \$2}')
        diaTracer: \$(java -jar ${tools_folder}/diaTracer-*.jar | head -n1 | awk -F- '{print \$2}')
    END_VERSIONS
    """

}
