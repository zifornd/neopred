include { ARCASHLA_EXTRACT               } from '../../../modules/nf-core/arcashla/extract'
include { ARCASHLA_GENOTYPE              } from '../../../modules/local/arcashla/genotype'

//Workflow

workflow HLA_TYPING {
    take:
    samplesheet               // channel: [ samplesheet ]
    bam_files                 // channel: [ val(meta), [ bam_files ] ]
    after_br                  // channel: [ batch removed genesymbols ]
    batch                     // batch variable
    design                    // design variable
    patient_id                // patient ID variable


    main:
    ch_versions = Channel.empty()
    ch_fastq = Channel.empty()
    // 
    // Module: Extracts reads mapped to chromosome 6 and any HLA decoys or chromosome 6 alternates
    //  

    ARCASHLA_EXTRACT (bam_files)
    ch_fastq = ch_versions.collect(ARCASHLA_EXTRACT.out.extracted_reads_fastq)
    ch_versions = ch_versions.mix(ARCASHLA_EXTRACT.out.versions)

    //
    // Module: Genotypes HLA alleles from extracted reads
    //

    ARCASHLA_GENOTYPE (ch_fastq)
    ch_versions = ch_versions.mix(ARCASHLA_GENOTYPE.out.versions)

    //
    // Module: Merges genotyping output for multiple samples into a single json file
    //

    ARCASHLA_MERGE (ARCASHLA_GENOTYPE.out.gt_json.map{it[1]}.collect())
    ch_versions = ch_versions.mix(ARCASHLA_MERGE.out.versions)

    //
    // Module: Changes alleles from its input form to a specified P-group or G-group nomenclature 
    //         or a specified 1, 2 or 3 fields in resolution
    //

    ARCASHLA_CONVERT (ARCASHLA_MERGE.out.merged_gt)
    ch_versions = ch_versions.mix(ARCASHLA_CONVERT.out.versions)

    //
    // Module: Generates a png file consisting the frequency of HLA as a plot
    //

    ARCASHLA_PLOT (ARCASHLA_CONVERT.out.gt_group,samplesheet,after_br)
    ch_versions = ch_versions.mix(ARCASHLA_PLOT.out.versions)


    emit:
    genotype   = ARCASHLA_GENOTYPE.out.genes_json   // channel: [meta, genes]
    hla_plot   = ARCASHLA_PLOT.out.hla_plot         // channel: [hla_frequency_plot]
    versions   = ch_versions                        // channel: [ versions.yml ]

}
