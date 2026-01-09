# Haplotree Parser for low-coverage ancient DNA data

Given an input BAM-file (e.g. Human extracted Reads from quicksand), this pipeline re-maps these sequences to the RSRS, extracts damaged sequences based on C-to-T damage and creates summary statistics for each Haplogroup-Node in the PhyloTree 17 release.

This is NOT a haplogroup caller! The output tables/trees can help to decide, which nodes are best presented in your data.

## Requirements

- singularity
- nextflow v22.10 (or larger)

## Usage

please prepare INPUT-DIR, a folder with BAM files containing human sequences. Then run the pipeline with:

```
NXF_VER=24.04.4 nextflow run merszym/haplotree --split INPUT-DIR
```

## Output

The pipeline produces three files for each input-file (in `out/06_haplogroups/`):
1. Summary statistics for each node in the tree (**Complete**)
2. Summary statistics for each node in the tree that is covered by sequences (**Most useful**)
3. Summary statistics for a single 'best' path (**provides a quick look, which nodes are dominant**)

The columns in the tables are

- **PhyloTree:** The Haplogroup Relationship in PhyloTree 17
- **NodeCoverage:** Accumulated values of the 'PositionCoverage' column in the full node subtree (sum of all children)
- **PositionCoverage:** For all Haplogroup-defining positions, count the ones that share the diagnostic allel and the once covered by at least 1 read (see 'ReadCoverage'). E.g. 2/2 = 2 positions out of 2 share the haplogroup-defining state
- **Penalty:** If a node is not supported, the penalty shows the number of nodes that are skipped before reaching a supported child-node in the same branch. The higher the penalty, the more a branch is supported only by single hits in the tips, rather than by following the tree top-down
- **ReadCoverage:** Shows for each _diagnostic position_ the number of sequences covering that position and the support for the diagnostic state in %.

## Example

An example for the 'best-path' file of data from an haplogroup individual (~3x genomic coverage in the deaminated sequences) can be found [here](assets/example_hapU.tsv)

See the high penalty-values for the U-subgroups, indicating a better support for the ancestral node (U)


## Workflow

1. Files are (re-)mapped to the RSRS (Reconstructed Sapiens Reference Sequence) with *BWA*
    - Saved to: `out/01_bwa/`
2. Files are filtered for mapping quality (25) and sorted using *samtools*
3. PCR duplicates are removed with *bam-rmdup*
    - Saved to `out/02_uniq`
4. Sequences are removed that overlap low-complexity poly-c stretches in the mtDNA genome (positions 303-315, 513-576, 3565-3576 and 16184-16193) with *bedtools intersect*
    - Saved to `out/03_bedfilter`
5. Variable positions in the alignment are detected (positions, in which the RSRS-reference base has less than 50% support)
    - Saved to `out/04_pielup`
6. Sequences are extracted that have a C-to-T substitution in the first or last three bases (unless this subsitution is one of the defined variable positions)
    - Saved to `05_deaminated`
7. Set the mapping quality score of the first and last three T-bases to 0 (masking)
8. Create a pileup-version of the files using `samtools`
    - Saved to `05_deaminated`
9. Walk through the PhyloTree-file and create thee different summary statistics for each node, based on the created pileup
    - Saved to `06_haplogroups`
