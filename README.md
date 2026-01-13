# ALPHABET
**AL**PHABET **p**rovides **ha**plogroup **b**ranch **e**s**t**imates

ALPHABET is a tool for the detection of mtDNA haplogroup-branches in low coverage or highly contaminated ancient DNA samples. It uses the PhyloTree 17 phylogeny to provide coverage-statistics for each haplogroup-node to find the haplogroup branch(es) with the highest support.

ALPHABET attempts to find the _branches_ with the highest support, it is _not_ a haplogroup caller!

## Requirements

- singularity
- nextflow v22.10 (or larger)

## Usage

The input for ALPHABET is a folder with BAM-files. These BAM-files should contain human mtDNA sequences (e.g. the Hominidae 'Extracted Reads' provided by the [quicksand](https://github.com/mpieva/quicksand) pipeline).

ALPHABET maps these sequences to the RSRS, extracts deaminated sequences (based on C-to-T damage) and creates summary statistics for each Haplogroup-Node in the PhyloTree 17 release (see workflow section below).

### quickstart

(Work in Progress)

Run alphabet:

```
nextflow run merszym/alphabet --split split
```

## Output

The pipeline produces three files for each input-file (in `out/06_haplogroups/`):
1. Summary statistics for each node in the tree
2. Filtered tree, showing only nodes with 70% branch-support or higher
3. Filtered tree, showing only nodes with 70% branch-support or higher AND all branches removed that have 3 consecutive gaps (**Most confident!**)

### Columns explained

- **Parent**: The parent-haplogroup (because the tree is not easy to parse)
- **PhyloTree:** The Haplogroup Relationship in PhyloTree 17
- **Penalty:** The shortest distance to the a child-node with support. 0 if the node is supported, -1 if no children has support. 
- **BranchSupport:** Accumulated PositionSupport in the branch 
- **BranchSupportPercent:** Accumulated PositionSupport in the branch (in %)
- **PositionSupport:** The number for haplogroup-defining positions that are covered by sequences and share the required state (see 'ReadCoverage')
- **SequenceSupport:** Coverage support for each _diagnostic position_. Shows the number of sequences covering that position and support the diagnostic state.

## Workflow

### 1. Mapping with BWA

Files are (re-)mapped to the RSRS (Reconstructed Sapiens Reference Sequence) with *BWA* and saved to `out/01_bwa/`

### 2. Filter Alignment
Files are filtered for minimum mapping quality (25) and minimum length (35). Then the alignment is sorted using *samtools*

### 3. Duplicate Removal

PCR duplicates are removed with *bam-rmdup* and saved to `out/02_uniq/`

### 4. Remove Poly-C Stretches

Sequences are removed from the alignment that overlap low-complexity poly-c stretches (positions 303-315, 513-576, 3565-3576 and 16184-16193) with *bedtools intersect*. These poly-c stretches can introduce face C-to-T deamination patterns, e.g.: ![](assets/img/poly_c_stretches.png).

The filtered alignment is saved to `out/03_bedfilter/`

### 5. Mask variable positions 

Variable positions can cause C-to-T differences in the first and last 3 bases because of haplogroup differences, rather than DNA damage. Example: ![](assets/img/variable_positions.png)

Positions in the alignment are masked if the majority of sequences (but at least 2 sequences) shows a differenct base than the reference.

The pileup and the extracted positions are saved to `out/04_pileup/`. The positions are only masked for the extraction of deaminated sequences.

### 6. Extract deaminated sequences

Sequences are extracted that have a C-to-T substitution in the first or last three bases (unless this subsitution is one of the masked positions). Deaminated sequences are saved to `05_deaminated/`

### 7. Mask first and last three bases 

Set the mapping quality score of the first and last three T-bases of all deaminated sequences to 0. They are ignored in the phylotree-analysis

### 8. Haplogroup Statistics

Walk through the PhyloTree-file and create summary statistics for each haplogroup node. The resulting tables are saved to   `06_haplogroups/`

# Resources

This repository includes the RSRS-based PhyloTree17 XML-file provided under MIT License by the [Institute of Genetic Epidemiology, Insbruck](https://github.com/genepi/phylotree-rsrs-17/blob/main/src/tree.xml)