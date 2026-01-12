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
2. Summary statistics for the tree, but branches with penalty-values >=3 or 0% node support are pruned
3. Summary statistics for a single 'best' path (**provides a quick look, which nodes are dominant**)
    - Note: This is experimental! This 'best' path is often too simple and needs manual verification!

The columns in the tables are

- **PhyloTree:** The Haplogroup Relationship in PhyloTree 17
- **NodeCoverage:** Accumulated values of the 'PositionCoverage' column in the full node subtree (sum of all children)
- **PositionCoverage:** For all Haplogroup-defining positions, count the ones that share the diagnostic allel and the once covered by at least 1 read (see 'ReadCoverage'). E.g. 2/2 = 2 positions out of 2 share the haplogroup-defining state
- **Penalty:** If a node is not supported, the penalty shows the closest distance to the a child-node with support in the same branch.
- **ReadCoverage:** Shows for each _diagnostic position_ the number of sequences covering that position and the support for the diagnostic state in %.

## Example

An example for the 'best-path' file of data from an haplogroup individual (~3x genomic coverage in the deaminated sequences) can be found [here](assets/example_hapU.tsv). The interpretation:

- High penalty-value (3) for the U8 subgroup indicates a better support for the ancestral node (U or U2'3'4'7'8'9)
- The U8 is only supported by hits to K-positions
- The detected 'K' and 'K1a4a1e' haplogroup positions are only supported by 'ancestral' positions (16311T! and 152!, see L2'3'4'5'6 and L3)


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
