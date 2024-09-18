process RSEQC_SIZE_DOWNSAMPLING {
    tag "$meta.id"
    label 'process_medium'

    //conda "${moduleDir}/environment.yml"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
    //    'biocontainers/rseqc:5.0.3--py39hf95cd2a_0' }"

    input:
    //tuple val(meta), path(bam), path(bai)
    //path  bed
    tuple val(meta), path(txt)

    output:
    //tuple val(meta), path("*.geneBodyCoverage.r"), emit: rscript
    //path("*.geneBodyCoverage.txt")               , emit: gbc_txt
    //path("*.geneBodyCoverage.curves.pdf")        , emit: pdf_curves
    tuple val(meta), path("*_stat_tmp.txt"), emit: stat_tmp
    path  "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //gene_body_cvg_qc:- gives RNA-seq read coverage over the gene body
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ds_check_size.sh \\
        $txt \\
        "${prefix}_stat_tmp.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | sed -e 's/, .*//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gene_body_coverage.r

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | sed -e 's/, .*//g')
    END_VERSIONS
    """
}

