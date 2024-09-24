process RSEQC_BAM_DOWNSAMPLING {
    tag "$meta.id"
    label 'process_medium'

    //conda "${moduleDir}/environment.yml"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
    //    'biocontainers/rseqc:5.0.3--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    //path  bed
    tuple val(meta), path(txt)

    output:
    //tuple val(meta), path("*.geneBodyCoverage.r"), emit: rscript
    //path("*.geneBodyCoverage.txt")               , emit: gbc_txt
    //path("*.geneBodyCoverage.curves.pdf")        , emit: pdf_curves
    tuple val(meta), path("*_downsampling.bam"), emit: downsample_bam
    path  "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //gene_body_cvg_qc:- gives RNA-seq read coverage over the gene body
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    COUNT="\$(ds_check_size.sh $txt)"
    downsampling.sh $bam "\${COUNT}"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}._downsampling.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
