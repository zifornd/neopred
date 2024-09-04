include { ARCASHLA_EXTRACT               } from '../../../modules/nf-core/arcashla/extract'
include { ARCASHLA_GENOTYPE              } from '../../../modules/local/arcashla/genotype'

//Workflow

workflow HLA_TYPING {
    take:
    bam_files                 // channel: [ val(meta), [ bam_files ] ]


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

    emit:

    genotype   = ARCASHLA_GENOTYPE.out.gt_json   // channel: [meta, genoty]
    versions   = ch_versions                     // channel: [ versions.yml ]

}
