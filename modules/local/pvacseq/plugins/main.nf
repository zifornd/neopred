process PVACTOOLS_PVACSEQPLUGINS {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/griffithlab/pvactools:4.1.1' :
        'griffithlab/pvactools:4.1.1' }"

    input:
    tuple val(meta), path(vcf)
    val algorithms

    output:
    tuple val(meta), path("*.tsv")  , optional:true, emit: results

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.filter"
    def callers = algorithms ?: "MHCflurry,MHCnuggetsII"
    """
    pvacseq run \\
        $vcf {params.tumor} {params.HLA} $callers {params.output_dir} \\
        -e1 {params.epitope_lengths} \\
        -t $task.cpus \\
        -r 5 --pass-only \\
        --fasta-size 200 \\
        -d 250 \\
        --iedb-install-directory /opt/iedb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pvacseq: \$(pvacseq --version | sed 's/Python //g')
    END_VERSIONS
    """
}