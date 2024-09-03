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

workflow PRE_VARIANTCALLING{

    take:
    bam
    bam_bai
    fasta
    fasta_fai
    known_indels      // channel: /path/to/reference.vcf.gz
    known_indels_tbi  // channel: /path/to/reference.vcf.gz.tbi

    //bam               // channel: [ val(meta), [ bamfiles ] ]
    //fasta             // channel: /path/to/genome.fa
    //fai               // channel: /path/to/genome.fai
    //dict              // channel: /path/to/genome.dict
    //markdup_bam       // channel: []
    //markdup_bai

    main:

    ch_versions = Channel.empty()

    PICARD_COLLECTMULTIPLEMETRICS(bam_bai,fasta,fasta_fai)
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())

    PICARD_ADDORREPLACEREADGROUPS (bam, [[:],[]], [[:],[]])
    ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions.first())

    //SAMTOOLS_FLAGSTAT_BAMRG (PICARD_ADDORREPLACEREADGROUPS.out.bam)

    PICARD_MARKDUPLICATES (PICARD_ADDORREPLACEREADGROUPS.out.bam, fasta,fasta_fai)
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    //SAMTOOLS_FLAGSTAT_BAMMARKDUB (PICARD_MARKDUPLICATES.out.bam)

    ch_bam_bai_int = PICARD_MARKDUPLICATES.out.bam
        .join(PICARD_MARKDUPLICATES.out.bai, by: [0])
        .map{ meta, bam, bai -> [ meta, bam, bai, [] ] }

    PICARD_CREATESEQUENCEDICTIONARY (fasta)
    ch_dict = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
    ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions.first())

    GATK4_SPLITNCIGARREADS (ch_bam_bai_int,fasta,fasta_fai,ch_dict)
    ch_versions = ch_versions.mix(GATK4_SPLITNCIGARREADS.out.versions.first())

    GATK4_BASERECALIBRATOR(ch_bam_bai_int, fasta.map{ meta, it -> it }, fasta_fai.map{ meta, it -> it }, ch_dict.map{ meta, it -> it }, known_indels, known_indels_tbi)
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first())

    ch_input_bqsr = ch_bam_bai_int.join(GATK4_BASERECALIBRATOR.out.table,by:[0])
        .map{ meta, bam, bai, empty, table -> [ meta, bam, bai, table, empty] }

    GATK4_APPLYBQSR(ch_input_bqsr, fasta.map{ meta, it -> it }, fasta_fai.map{ meta, it -> it }, ch_dict.map{ meta, it -> it })
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions.first())

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
