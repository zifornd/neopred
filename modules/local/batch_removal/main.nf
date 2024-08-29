process BATCH_REMOVAL {
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-limma_bioconductor-sva_r-base_r-optparse:d938a560bf123c62' :
        'community.wave.seqera.io/library/bioconductor-limma_bioconductor-sva_r-base_r-optparse:d938a560bf123c62' }"


    input:
    path (samplesheet)
    path (tpm_gene)

    output:
    path "*_tpm.genesymbol.batchremoved.csv", emit: after_br
    path "*_tpm.genesymbol.csv"             , emit: before_br
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/

    """
    batch_removal.R \\
        -e ${tpm_gene} \\
        -c batch \\
        -d design \\
        -m ${samplesheet} \\
        -b design_batch_tpm.genesymbol.csv \\
        -a design_batch_tpm.genesymbol.batchremoved.csv
    cp design_batch_tpm.genesymbol.batchremoved.csv tpm.genesymbol.batchremoved.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
