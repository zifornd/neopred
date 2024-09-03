//
// Pre-processing Base Quality Recalibration with GATK4
//

include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS } from '../../../modules/local/picard/collectmetrics/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_ADDORREPLACEREADGROUPS } from '../../../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_BAMRG     } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_BAMMARKDUB } from '../../../modules/nf-core/samtools/flagstat/main'
include { PICARD_CREATESEQUENCEDICTIONARY               } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { GATK4_SPLITNCIGARREADS                        } from '../../../modules/nf-core/gatk4/splitncigarreads/main'
include { GATK4_BASERECALIBRATOR } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR              } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { GATK4_INDEXFEATUREFILE } from '../../../modules/nf-core/gatk4/indexfeaturefile/main'
/*
include { PICARD_CREATESEQUENCEDICTIONARY               } from '../../modules/nf-core/picard/createsequencedictionary/main'
include { GATK4_SPLITNCIGARREADS                        } from '../../modules/nf-core/gatk4/splitncigarreads/main'
include { GATK4_INDEXFEATUREFILE                     } from '../../modules/nf-core/gatk4/indexfeaturefile/main'
include { GATK4_BASERECALIBRATOR } from '../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR              } from '../../modules/nf-core/gatk4/applybqsr/main'
*/
workflow PRE_VARIANTCALLING{

    take:

    bam
    fasta
    fai
    intervals //channel: [mandatory] [ intervals, num_intervals ] (or [ [], 0 ] if no intervals)

    //bam               // channel: [ val(meta), [ bamfiles ] ]
    //fasta             // channel: /path/to/genome.fa
    //fai               // channel: /path/to/genome.fai
    //dict              // channel: /path/to/genome.dict
    //markdup_bam       // channel: []
    //markdup_bai
    //known_indels      // channel: /path/to/reference.vcf.gz
    //known_indels_tbi  // channel: /path/to/reference.vcf.gz.tbi

    main:

    ch_versions = Channel.empty()

    //
    //PICARD_COLLECTMULTIPLEMETRICS
    //
    PICARD_COLLECTMULTIPLEMETRICS(bam,fasta,fai)

    PICARD_ADDORREPLACEREADGROUPS (bam, [[:],[]], [[:],[]])

    //SAMTOOLS_FLAGSTAT_BAMRG (PICARD_ADDORREPLACEREADGROUPS.out.bam)

    PICARD_MARKDUPLICATES (PICARD_ADDORREPLACEREADGROUPS.out.bam, fasta,fai)

    //SAMTOOLS_FLAGSTAT_BAMMARKDUB (PICARD_MARKDUPLICATES.out.bam)

    PICARD_CREATESEQUENCEDICTIONARY (fasta.join(fai,by:[0]))

    ch_dict = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
    ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)

    ch_bam_bai = PICARD_MARKDUPLICATES.out.bam
        .join(PICARD_MARKDUPLICATES.out.bai, by: [0])
        .join(intervals, by: [0])
        .map {
        meta, bam, bai, intervals -> {[meta, bam, bai, intervals]}
        }
        .set{ch_bam_bai_int}

    ch_bam_bai_int.view()


    GATK4_SPLITNCIGARREADS (PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai,by:[0]),fasta,fai,ch_dict)
    //
    // Creates index for given vcf file. - Can be used when .tbi file of known indels are not provided
    //

    //
    // Create sequence dictionary file for reference fasta
    //
    /*
    PICARD_CREATESEQUENCEDICTIONARY (
        fasta
    )

    ch_dict = PICARD_CREATESEQUENCEDICTIONARY.out.dict
    ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)

    markdup_bam
    .branch {
            meta, files ->
                wes: meta.seq_type == 'WES'
                    return [ meta, files ]
                rna: meta.seq_type == 'RNA'
                    return [ meta, files ]
    }
    .set { markdup_bam }

    markdup_bai
    .branch {
            meta, files ->
                wes: meta.seq_type == 'WES'
                    return [ meta, files ]
                rna: meta.seq_type == 'RNA'
                    return [ meta, files ]
    }
    .set { markdup_bai }

    //
    // MODULE: Splits reads that contain Ns in their cigar string
    //
    GATK4_SPLITNCIGARREADS (
        markdup_bam.rna,
        fasta,
        fai,
        ch_dict
        )
        .bam
        .mix(markdup_bam.wes)
        .set { ch_splitreads_bam }

    ch_splitreads_bai = GATK4_SPLITNCIGARREADS.out.bai.mix(markdup_bai.wes)

    //ch_splitreads = GATK4_SPLITNCIGARREADS.out.bam

    //
    // Generate recalibration table for BQSR
    //
    GATK4_BASERECALIBRATOR (ch_splitreads_bam, fasta, fai, ch_dict, known_indels, known_indels_tbi )
    ch_bqsr_table    = GATK4_BASERECALIBRATOR.out.table
    ch_bqsr_versions = GATK4_BASERECALIBRATOR.out.versions

    //
    // Apply BQSR to BAM file
    //
    GATK4_APPLYBQSR ( ch_splitreads_bam, ch_bqsr_table, fasta, fai, ch_dict)
    ch_bqsr_bam   = GATK4_APPLYBQSR.out.bam
    ch_bqsr_bai   = GATK4_APPLYBQSR.out.bai
    ch_bqsr_cram  = GATK4_APPLYBQSR.out.cram
    ch_versions   = ch_versions.mix(GATK4_SPLITNCIGARREADS.out.versions,GATK4_BASERECALIBRATOR.out.versions,GATK4_APPLYBQSR.out.versions)
    */

    emit:
    /*
    bqsr_table     = ch_bqsr_table              // channel: [ val(meta), bqsr_table      ]
    bqsr_bam       = ch_bqsr_bam                // channel: [ val(meta), bqsr_bam        ]
    bqsr_bai       = ch_bqsr_bai                // channel: [ val(meta), bqsr_bai        ]
    bqsr_cram      = ch_bqsr_cram               // channel: [ val(meta), bqsr_cram       ]
    splitreads_bam = ch_splitreads_bam         // channel: [ val(meta), splitreads      ]
    splitreads_bai = ch_splitreads_bai         // channel: [ val(meta), splitreads      ]
    dict	   = ch_dict		       // channel: [ val(meta), dict		]
    */
    versions       = ch_versions                // channel: [ versions.yml		]
}
