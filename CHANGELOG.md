## v0.8 [02.02.2026]
### Changes
- Add a first version of Sima de los huesos and Denisovan diagnostic positions to the tree XML
- Add a upper hierarchie (NA'SIMA'DEN) for the common postions of Neanderthal, Sima and Denisova 
- Add a first version of a final report
- Add stats and output-files for Deduped and Deaminated sequences separately
- Remove the `min_support` and `max_gaps` flags, as they are not helpful in mixed samples (that was the original intention...)

## v0.7 [30.01.2026]
### Changes
- add Neanderthal haplogroups to the phylotree XML file, as published in [Andreeva et al. 2022](https://www.nature.com/articles/s41598-022-16164-9)
- Update LICENSE to include the MIT LICENSE of the [Institute of Genetic Epidemiology, Insbruck](https://github.com/genepi/phylotree-rsrs-17/)
- Include the `--penalty_plus` flag (default: 2) to display nodes with higher penalty than the lowest.

## v0.6 [19.01.2026]
### Changes
- add `TotalMismatch` column, indicating the total number of mismatches in the branch support between covered and supported positions
- add `SumOfGaps` column, that shows the accumulated number of gaps in this branch
- add `DistanceToBest` column that shows the difference of supported branch positions to the maximum number of supported positions in any branch
- update the `Penalty` column to show `SumOfGaps` + `TotalMismatch` + `DistanceToBest`
- Rename (and add) output tree-tables
    - rename the full table from `all_groups` to `raw` 
    - add a table that shows the path to the nodes with the best two (minimum) penalty scores (`best`)

## v0.5 [14.01.2026]
### Changes
- Add the `RequiredGaps` column in the summary that, looking back, counts how many intermediate nodes were skipped to get there
    (different to the `Penalty` column that looks forward to the closest child)
- introduce the `--max_gaps` and `--min_support` flags to filter the tree based on the BranchSupport and the RequiredGaps columns


## v0.4 [12.01.2026]
### Changes
- Reformat parse_pylotree code (more readable and maintainable)
- Use BranchSupport as a new filter-metric
- Update the output-files (All Nodes, 70% Support, 70% Support + Max 3 Gaps)
- Add MIT LICENSE

## v0.3
### Changes
- Call a position with ! (re-mutation) only as 'found', if the same position was not requested before. This is necessary, because they hit tip-positions way too often and accumulate hits in large branches. 
    - e.g. 16311T in L3 -> 16311T! in U1b3 (not found)
    - e.g. 16311T in L3 -> 16311A! in U5? (found)
- Second output (groups with coverage) now clips all nodes that have a penalty >= 3
    - this _can_ make the second file less useful in some instances, but much more useful in many more :)