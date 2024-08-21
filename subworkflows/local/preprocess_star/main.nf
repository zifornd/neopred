//
// Alignment with STAR and getting the Stats
//

include { STAR_ALIGN              } from '../../../modules/nf-core/star/align'
include { SAMTOOLS_SORT           } from '../../../modules/nf-core/samtools/sort'
include { SAMTOOLS_INDEX          } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_STATS          } from '../../../modules/nf-core/samtools/stats'


workflow PREPROCESS_STAR {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    index               // channel: [ val(meta), [ index ] ]
    gtf                 // channel: [ val(meta), [ gtf ] ]
    star_ignore_sjdbgtf // boolean: when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file 
    seq_platform        // string : sequencing platform
    seq_center          // string : sequencing center
    is_aws_igenome      // boolean: whether the genome files are from AWS iGenomes
    fasta               // channel: /path/to/fasta

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with STAR
    //
    ch_orig_bam       = Channel.empty()
    ch_log_final      = Channel.empty()
    ch_log_out        = Channel.empty()
    ch_log_progress   = Channel.empty()
    ch_bam_sorted     = Channel.empty()
    ch_bam_transcript = Channel.empty()
    ch_fastq          = Channel.empty()
    ch_tab            = Channel.empty()

    //
    // Module: Star Align
    //
    STAR_ALIGN ( reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
    ch_orig_bam       = STAR_ALIGN.out.bam
    ch_log_final      = STAR_ALIGN.out.log_final
    ch_log_out        = STAR_ALIGN.out.log_out
    ch_log_progress   = STAR_ALIGN.out.log_progress
    ch_bam_sorted     = STAR_ALIGN.out.bam_sorted
    ch_bam_transcript = STAR_ALIGN.out.bam_transcript
    ch_fastq          = STAR_ALIGN.out.fastq
    ch_tab            = STAR_ALIGN.out.tab
    ch_versions       = ch_versions.mix(STAR_ALIGN.out.versions.first())

    //
    // Module: Samtools sort
    //
    SAMTOOLS_SORT ( ch_orig_bam, fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    //
    // Module: Samtools index
    //
    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    //
    // Module: Samtools stats
    //
    SAMTOOLS_SORT.out.bam
    .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
    .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
    .map {
        meta, bam, bai, csi ->
            if (bai) {
                [ meta, bam, bai ]
            } else {
                [ meta, bam, csi ]
            }
    }
    .set { ch_bam_bai }
    SAMTOOLS_STATS ( ch_bam_bai, fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)  

    emit:
    orig_bam       = ch_orig_bam                    // channel: [ val(meta), bam            ]
    log_final      = ch_log_final                   // channel: [ val(meta), log_final      ]
    log_out        = ch_log_out                     // channel: [ val(meta), log_out        ]
    log_progress   = ch_log_progress                // channel: [ val(meta), log_progress   ]
    bam_sorted     = ch_bam_sorted                  // channel: [ val(meta), bam_sorted     ]
    bam_transcript = ch_bam_transcript              // channel: [ val(meta), bam_transcript ]
    fastq          = ch_fastq                       // channel: [ val(meta), fastq          ]
    tab            = ch_tab                         // channel: [ val(meta), tab            ]

    bam_sort    = SAMTOOLS_SORT.out.bam          // channel: [ val(meta), [ bam ] ]
	 stats          = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), [ stats ] ]

    versions       = ch_versions                    // channel: [ versions.yml ]
}

