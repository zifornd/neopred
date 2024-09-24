process ARCASHLA_PLOT {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-complexheatmap_r-dplyr_r-ggplot2_r-optparse_pruned:3526b6f9e2287908' :
        'community.wave.seqera.io/library/bioconductor-complexheatmap_r-dplyr_r-ggplot2_r-optparse_pruned:3526b6f9e2287908' }"


    input:
    path gt
    path samplesheet
    path after_br
    val batch
    val design
    val patid

    output:
    path ("*.pdf")       , emit: hla_plot
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = "0.5.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    hla_plot.R \\
        --hla $gt \\
        --meta $samplesheet \\
        --expression $after_br \\
        --design $design \\
        --patID $patid \\
        --source ${baseDir}/bin/hlaonco_plot.R \\
        --outdir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
