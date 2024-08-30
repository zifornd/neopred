//
// Quantification with Salmon
//

include { SALMON_QUANT   } from '../../../modules/nf-core/salmon/quant'
include { TX2GENE        } from '../../../modules/local/tx2gene'
include { TXIMPORT       } from '../../../modules/local/tximport'

workflow QUANTIFY_SALMON {
    take:
    bam_files                 // channel: [ val(meta), [ bam_files ] ]
    index                     // channel: /path/to/salmon/index/
    transcript_fasta          // channel: /path/to/transcript.fasta
    gtf                       // channel: /path/to/genome.gtf
    alignment_mode            //    bool: Run Salmon in alignment mode
    lib_type                  //     val: String to override Salmon library type

    main:

    ch_versions = Channel.empty()

    //
    // Quantify and merge counts across samples
    //
    // NOTE: MultiQC needs Salmon outputs
    SALMON_QUANT ( bam_files, index, transcript_fasta, alignment_mode, lib_type )
    ch_results = SALMON_QUANT.out.results
    ch_multiqc = ch_results
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    TX2GENE ( ch_results.collect{it[1]}, 'salmon', gtf )
    ch_versions = ch_versions.mix(TX2GENE.out.versions)

    TXIMPORT ( ch_results.collect{it[1]}, TX2GENE.out.tsv.collect(), 'salmon' )
    ch_versions = ch_versions.mix(TXIMPORT.out.versions)

    emit:
    results                       = ch_results                      // channel: [ val(meta), results_dir ]
    multiqc                       = ch_multiqc                      // channel: [ val(meta), files_for_multiqc ]

    tpm_gene                      = TXIMPORT.out.tpm_gene           //    path *gene_tpm.tsv

    versions                      = ch_versions                     // channel: [ versions.yml ]
}
