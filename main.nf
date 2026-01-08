// Import modules

include { MAP_BWA            } from './modules/local/map_bwa'
include { BEDTOOLS_INTERSECT } from './modules/local/bedtools_intersect'
include { SAMTOOLS_FQ2BAM    } from './modules/local/samtools_fq2bam'
include { SAMTOOLS_SORT      } from './modules/local/samtools_sort'
include { BAM_RMDUP          } from './modules/local/bam_rmdup'
include { SAMTOOLS_MPILEUP   } from './modules/local/samtools_mpileup'
include { SAMTOOLS_MPILEUP_DEAM3  } from './modules/local/samtools_mpileup'
include { GET_VARIABLE_POSITIONS  } from './modules/local/get_variable_positions'
include { FILTER_BAM         } from './modules/local/filter_bam'
include { MASK_DEAMINATION   } from './modules/local/mask_deamination'

// load the files

ch_split      = Channel.fromPath("${params.split}/*"  ,checkIfExists:true) // input-data
ch_reference  = Channel.fromPath("${params.reference}" ,checkIfExists:true) // Reference genome (RSRS)
ch_bedfile    = Channel.fromPath("${params.bedfile}"  ,checkIfExists:true) // bedfile (poly-c stretches)

ch_versions = Channel.empty()

// some required functions
def has_ending(file, extension){
    return extension.any{ file.toString().toLowerCase().endsWith(it) }
}

//
// MAIN WORKFLOW
//

workflow {

// add a first meta
ch_split.map{it -> [['sample': it.baseName, 'id':it.baseName], it] }.set{ ch_split }

//split input into bam- and fastq-files
ch_split.branch {
    bam: it[1].getExtension() == 'bam' 
    fastq: has_ending( it[1], ["fastq","fastq.gz","fq","fq.gz"]) //input FASTQ need to be converted to bam-files
    fail: true
}
.set{ ch_split }

//
// 0. Fastq to BAM
//

SAMTOOLS_FQ2BAM(ch_split.fastq)

ch_versions = ch_versions.mix(SAMTOOLS_FQ2BAM.out.versions.first())

ch_split_bam = ch_split.bam.mix(
    SAMTOOLS_FQ2BAM.out.bam
)


//
// 1. Map BWA and sort
//

ch_split_for_bwa = ch_split_bam.combine(ch_reference)

MAP_BWA(ch_split_for_bwa)
ch_versions = ch_versions.mix(MAP_BWA.out.versions.first())


SAMTOOLS_SORT(MAP_BWA.out.bam)
ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

//
// 2. Deduplicate 
//

BAM_RMDUP( SAMTOOLS_SORT.out.bam )

// Add the Number of Unique Sequences to the meta
// first parse the bam-rmdup output file
stats = BAM_RMDUP.out.txt.splitCsv(sep:'\t', header:true, limit:1)

//now add the stats to the meta
BAM_RMDUP.out.bam
.combine(stats, by: 0)
.map{ meta, bam, stats ->
    [
        meta+["ReadsDeduped":stats["out"].replace(",","") as int],
        bam
    ]
}
.set{ ch_split_for_bed }
ch_versions = ch_versions.mix(BAM_RMDUP.out.versions.first())

//
// 3. Remove poly-c sequences
//

ch_split_for_bed = ch_split_for_bed.combine(ch_bedfile)

BEDTOOLS_INTERSECT(ch_split_for_bed)
ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions.first())


//
// 4. Make pileup
//

ch_pileup = BEDTOOLS_INTERSECT.out.bam.combine(ch_reference)

SAMTOOLS_MPILEUP(ch_pileup)

tsv = SAMTOOLS_MPILEUP.out.tsv

//
// 5. Get variable positions
//

GET_VARIABLE_POSITIONS(tsv)
ch_versions = ch_versions.mix(GET_VARIABLE_POSITIONS.out.versions.first())

//
// 6.Extract Deaminated Sequences
//

ch_for_deam3 = BEDTOOLS_INTERSECT.out.bam.combine(GET_VARIABLE_POSITIONS.out.tsv, by:0)

FILTER_BAM(ch_for_deam3)

//
// 7.Mask Deamination for the pileup
//

MASK_DEAMINATION(FILTER_BAM.out.bam)

//
// 8. Pileup the deaminated reads
//

SAMTOOLS_MPILEUP_DEAM3(MASK_DEAMINATION.out.bam)

//
// 
//

}