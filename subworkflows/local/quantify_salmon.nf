//
// Quantification with Salmon
//

include { SALMON_QUANT   } from '../../../modules/nf-core/salmon/quant'
include { TX2GENE        } from '../../../modules/local/tx2gene'
include { TXIMPORT       } from '../../../modules/local/tximport'

workflow QUANTIFY_SALMON {
    take:
    reads                     // channel: [ val(meta), [ reads ] ]
    index                     // channel: /path/to//index/
    transcript_fasta          // channel: /path/to/transcript.fasta
    gtf                       // channel: /path/to/genome.gtf
    lib_type                  //     val: String to override Salmon library type

    main:

    ch_versions = Channel.empty()

    //
    // Quantify and merge counts across samples
    //
    // NOTE: MultiQC needs Salmon outputs
    SALMON_QUANT ( reads, index, gtf, transcript_fasta, alignment_mode, lib_type )
    ch_results = SALMON_QUANT.out.results
    ch_multiqc = ch_results
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    TX2GENE ( ch_results.collect{it[1]}, pseudo_aligner, gtf )
    ch_versions = ch_versions.mix(TX2GENE.out.versions)

    TXIMPORT ( ch_results.collect{it[1]}, TX2GENE.out.tsv.collect(), pseudo_aligner )
    ch_versions = ch_versions.mix(TXIMPORT.out.versions)

    emit:
    results                       = ch_results                      // channel: [ val(meta), results_dir ]
    multiqc                       = ch_multiqc                      // channel: [ val(meta), files_for_multiqc ]

    tpm_gene                      = TXIMPORT.out.tpm_gene                  //    path *gene_tpm.tsv
    counts_gene                   = TXIMPORT.out.counts_gene               //    path *gene_counts.tsv
    lengths_gene                  = TXIMPORT.out.lengths_gene              //    path *gene_lengths.tsv
    counts_gene_length_scaled     = TXIMPORT.out.counts_gene_length_scaled //    path *gene_counts_length_scaled.tsv
    counts_gene_scaled            = TXIMPORT.out.counts_gene_scaled        //    path *gene_counts_scaled.tsv
    tpm_transcript                = TXIMPORT.out.tpm_transcript            //    path *gene_tpm.tsv
    counts_transcript             = TXIMPORT.out.counts_transcript         //    path *transcript_counts.tsv
    lengths_transcript            = TXIMPORT.out.lengths_transcript        //    path *transcript_lengths.tsv

    merged_counts_transcript      = TXIMPORT.out.counts_transcript         //    path: *.transcript_counts.tsv
    merged_tpm_transcript         = TXIMPORT.out.tpm_transcript            //    path: *.transcript_tpm.tsv

    versions                      = ch_versions                            // channel: [ versions.yml ]
}
