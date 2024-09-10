process PVACTOOLS {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/griffithlab/pvactools:4.1.1' :
        'griffithlab/pvactools:4.1.1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.tsv")  , optional:true, emit: results

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.filter"
    """
    pvacseq run \\
        {input.vcf} {params.tumor} {params.HLA} {params.callers} {params.output_dir} \\
        -e1 {params.epitope1_lengths} \\
        -t {threads} \\
        -r 5 --pass-only \\
        --fasta-size 200 \\
        -d 250 \\
        --iedb-install-directory {params.iedb}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
