process PCA_SAMPLE_CLUSTERING {
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-ggfortify_r-ggplot2_r-ggpubr_pruned:c2fbe12fb4826e03' :
        'community.wave.seqera.io/library/r-base_r-ggfortify_r-ggplot2_r-ggpubr_pruned:c2fbe12fb4826e03' }"


    input:
    path samplesheet
    val batch
    val design
    path after_br
    path before_br

    output:
    path "*_pca_plot_after.pdf" , emit: after_pca
    path "*_pca_plot_before.pdf", emit: before_pca
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pca.R \\
        -b ${before_br} \\
        -a ${after_br} \\
        -m ${samplesheet} \\
        -c ${batch} \\
        -g ${design} \\
        -i ${design}_${batch}_pca_plot_before.pdf \\
        -j ${design}_${batch}_pca_plot_after.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
