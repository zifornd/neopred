process PVACTOOLS_PVACSEQ {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/griffithlab/pvactools:4.1.1' :
        'griffithlab/pvactools:4.1.1' }"

    input:
    tuple val(meta), path(vcf)
    path(hla)
    val algorithms
    val neoantigen_epitope1_lengths

    output:
    tuple val(meta), path("*.tsv")  , optional:true, emit: results

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.filter"
    def callers = algorithms ?: "MHCflurry,MHCnuggetsII"
    """
    # Input arguments
    sample_name=$prefix
    # file="data.txt"  # Replace with your actual file path

    # Check if the sample exists in the file and dynamically extract all columns except the first one
    values=$(awk -v sample="$sample_name" '$1 == sample {for(i=2; i<=NF; i++) printf "%s%s", $i, (i<NF ? "," : "\n")}' "${hla}")

    # If no values were found, print null, otherwise print the values
    if [[ -z "$values" ]]; then
    echo "null"
    else
    echo "$values"
    pvacseq run \\
        $vcf ${meta.id} $values $callers ./ \\
        -e1 ${neoantigen_epitope1_lengths} \\
        -t $task.cpus \\
        -r 5 --pass-only \\
        --fasta-size 200 \\
        -d 250 \\
        --iedb-install-directory /opt/iedb
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pvacseq: \$(pvacseq --version | sed 's/Python //g')
    END_VERSIONS
    """
}
