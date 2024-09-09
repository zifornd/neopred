//Runs batch removal using limma method and generates PCA plot for clustering samples

include { BATCH_REMOVAL              } from '../../../modules/local/batch_removal'
include { PCA_SAMPLE_CLUSTERING      } from '../../../modules/local/pca_sample_clustering'

//workflow
workflow BATCH_REMOVAL_ANALYSIS {
    take:
    sample               // channel: [ val(meta), [ reads ] ]
    batch                // batch variable
    design               // design variable
    tpm_gene             // channel: /path/to/tx2gene/tpm_gene

    main:

    ch_versions = Channel.empty()

    //
    // Module: Batch removal
    //
    BATCH_REMOVAL ( sample,batch,design,tpm_gene )
    ch_versions       = ch_versions.mix(BATCH_REMOVAL.out.versions)

    PCA_SAMPLE_CLUSTERING ( sample,batch,design,BATCH_REMOVAL.out.after_br, BATCH_REMOVAL.out.before_br )
    ch_versions       = ch_versions.mix(PCA_SAMPLE_CLUSTERING.out.versions)

    emit:
    before_br_pca = PCA_SAMPLE_CLUSTERING.out.before_pca  // channel: [ pca_before_batch_removal_plot ]
    after_br_pca  = PCA_SAMPLE_CLUSTERING.out.after_pca   // channel: [ pca_after_batch_removal_plot ]

    versions      = ch_versions                           // channel: [ versions.yml ]
}
