// Import modules

include { MAP_BWA            } from './modules/local/map_bwa'
include { BEDTOOLS_INTERSECT } from './modules/local/bedtools_intersect'
include { SAMTOOLS_FQ2BAM    } from './modules/local/samtools_fq2bam'
include { SAMTOOLS_SORT      } from './modules/local/samtools_sort'
include { SAMTOOLS_COUNT     } from './modules/local/samtools_count'
include { BAM_RMDUP          } from './modules/local/bam_rmdup'
include { SAMTOOLS_MPILEUP   } from './modules/local/samtools_mpileup'
include { SAMTOOLS_MPILEUP_DEAM3  } from './modules/local/samtools_mpileup'
include { GET_VARIABLE_POSITIONS  } from './modules/local/get_variable_positions'
include { FILTER_BAM         } from './modules/local/filter_bam'
include { MASK_DEAMINATION   } from './modules/local/mask_deamination'
include { SUMMARIZE_PHYLOTREE as SUMMARIZE_PHYLOTREE_DEAM   } from './modules/local/parse_phylotree'
include { SUMMARIZE_PHYLOTREE as SUMMARIZE_PHYLOTREE_UNIQUE } from './modules/local/parse_phylotree'

// load the files

ch_split      = Channel.fromPath("${params.split}/*"   ,checkIfExists:true) // input-data
ch_reference  = Channel.fromPath("${params.reference}" ,checkIfExists:true) // Reference genome (RSRS)
ch_bedfile    = Channel.fromPath("${params.bedfile}"   ,checkIfExists:true) // bedfile (poly-c stretches)
ch_treexml    = Channel.fromPath("${params.treexml}"   ,checkIfExists:true) // phylotree17 XML

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
//include placeholder values for later overwriting
BAM_RMDUP.out.bam
.combine(stats, by: 0)
.map{ meta, bam, stats ->
    [
        meta+[
            "Sequences":"Deduped",
            "ReadsMapped":stats["in"].replace(",","") as int, 
            "ReadsDeduped":stats["out"].replace(",","") as int,
            "ReadsBedfiltered":"-",
            "ReadsDeam":"-"
            ],
        bam
    ]
}
.set{ ch_unique_bam }

ch_versions = ch_versions.mix(BAM_RMDUP.out.versions.first())

//
// 3. Remove poly-c sequences
//

ch_unique_bam = ch_unique_bam.combine(ch_bedfile)

BEDTOOLS_INTERSECT(ch_unique_bam)
ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions.first())

// get the counts
SAMTOOLS_COUNT( BEDTOOLS_INTERSECT.out.bam )
// add the counts to the meta

ch_bedfiltered = SAMTOOLS_COUNT.out.bam.map { meta, bam, count ->
    [ meta+[ 'ReadsBedfiltered': count as int ], bam ]
}

//
// 4. Make pileup
//

ch_pileup = ch_bedfiltered.combine(ch_reference)

SAMTOOLS_MPILEUP(ch_pileup)

tsv = SAMTOOLS_MPILEUP.out.tsv

//
// 4.1 and create the stats for all sequences
//

ch_final_unique = Channel.empty()

if(params.include_deduped){
    ch_for_unique_phylotree = SAMTOOLS_MPILEUP.out.tsv.combine(ch_treexml)

    SUMMARIZE_PHYLOTREE_UNIQUE(ch_for_unique_phylotree)

    SUMMARIZE_PHYLOTREE_UNIQUE.out.stats
        .map{ meta, txt ->
            def stats = txt.splitCsv(sep:'\t', header:true, limit:1)[0]
            meta+stats
        }
        .set{ ch_final_unique }
}

//
// 5. Get variable positions
//

GET_VARIABLE_POSITIONS(tsv)
ch_versions = ch_versions.mix(GET_VARIABLE_POSITIONS.out.versions.first())

//
// 6.Extract Deaminated Sequences
//

ch_for_deam3 = ch_bedfiltered.combine(GET_VARIABLE_POSITIONS.out.tsv, by:0)

FILTER_BAM(ch_for_deam3)

//
// 7.Mask Deamination for the pileup
//

MASK_DEAMINATION(FILTER_BAM.out.bam)

MASK_DEAMINATION.out.bam.map{
    meta, bam, count ->
    def deam = count.text.trim() as int
    [
        meta+['ReadsDeam':deam, "Sequences":"Deam"],
        bam
    ]
}.set{
    ch_deaminated
}

//
// 8. Pileup the deaminated reads
//

SAMTOOLS_MPILEUP_DEAM3(ch_deaminated)

//
// 9. Summarize the data based on the phylotree XML 
//

ch_for_phylotree = SAMTOOLS_MPILEUP_DEAM3.out.tsv.combine(ch_treexml)

SUMMARIZE_PHYLOTREE_DEAM(ch_for_phylotree)

//Add the node-stats to the meta
SUMMARIZE_PHYLOTREE_DEAM.out.stats
    .map{ meta, txt ->
        def stats = txt.splitCsv(sep:'\t', header:true, limit:1)[0]
        meta+stats
    }
    .set{ ch_final }


ch_final = ch_final.mix(ch_final_unique)

//
// 
// 10. Summary Report
//
//

header_map = [
'base' : ['File', 'Sequences', 'ReadsMapped','ReadsDeduped', 'ReadsBedfiltered', 'ReadsDeam'].join('\t'),
'hap' : ["Phylotree", "BranchSupport", "Penalty", "SumOfGaps", "SequenceSupport", "Note" ].join('\t')
]

//
// if the keys in the meta dont match the desired columns, map here the meta keys to the values...
//
value_map = [
    'base' : ['id', 'Sequences', 'ReadsMapped','ReadsDeduped', 'ReadsBedfiltered', 'ReadsDeam'].join('\t'),
]

def getVals = {String key, meta, res=[] ->
    if(value_map[key]) {
        header = value_map[key]
    } else {
        header = header_map[key]
    }
    // then
    header.split('\t').each{
        def entry_key = it.trim()
        res << meta[entry_key]
        }
    res.join('\t')
}

// Save the output
ch_final
    .collectFile( name:"final_report.tsv",
        seed:[
            header_map['base'],
            header_map['hap'],
        ].join('\t'), storeDir:".", newLine:true, sort:true
    ){
        [
            getVals('base', it),
            getVals('hap', it),
        ].join('\t')
    }
    .subscribe {
        println "[alphabet]: Summary reports saved"
    }

}