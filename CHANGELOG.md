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