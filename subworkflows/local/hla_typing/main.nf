// HLA Typing: Identify potential neoantigens

include { ARCASHLA_EXTRACT               } from '../../../modules/nf-core/arcashla/extract'
include { ARCASHLA_GENOTYPE              } from '../../../modules/local/arcashla/genotype'
include { ARCASHLA_MERGE                 } from '../../../modules/local/arcashla/merge'
include { ARCASHLA_CONVERT               } from '../../../modules/local/arcashla/convert'
include { ARCASHLA_PLOT                  } from '../../../modules/local/arcashla/plot'

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
    ch_fastq = ARCASHLA_EXTRACT.out.extracted_reads_fastq
    ch_versions = ch_versions.mix(ARCASHLA_EXTRACT.out.versions)

    //
    // Module: Genotypes HLA alleles from extracted reads
    //

    ARCASHLA_GENOTYPE (ch_fastq)
    ch_genotype = ARCASHLA_GENOTYPE.out.gt_json
    ch_versions = ch_versions.mix(ARCASHLA_GENOTYPE.out.versions)

    //
    // Module: Merges genotyping output for multiple samples into a single json file
    //

    ARCASHLA_MERGE (ch_genotype.collect{it[1]})
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
    ch_sample   = Channel.value(file(samplesheet)) 
    ARCASHLA_PLOT (ARCASHLA_CONVERT.out.gt_group,ch_sample,after_br,batch,design,patient_id)
    ch_versions = ch_versions.mix(ARCASHLA_PLOT.out.versions)


    emit:
    hla_log    = ARCASHLA_GENOTYPE.out.gt_log   // channel: [meta, genotype log]
    hla_plot   = ARCASHLA_PLOT.out.hla_plot         // channel: [hla_frequency_plot]
    versions   = ch_versions                        // channel: [ versions.yml ]

}
