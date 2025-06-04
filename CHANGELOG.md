## [0.1.1] - 2025-06-04

### Changed
- Updated the installation method.
- Modified the output directory structure.
- Added `-m` option to the `phylo` command.
- Fixed usage of `-eject` option in the `edid` command.
- Added `--species` option to the `edid` command.

### Added
- Added ploidy mode (`-w`) to the following commands: `pipeline`, `prepare`, `curate`, `phylo`, `marker`, and `primer`.
- In `prepare`, polyploid genomes now use `threads.csv` to filter input files.
- In `curate`, either `threads.tsv` or `threads.csv` is read to extract genes that match specified copy number patterns.
- In `phylo`, `polyphest` is used to generate `multree` and network (extended Newick format) as `species_tree`.
- In `marker`, `multree` is used as the `species_tree` for comparison with gene trees, and nRF values are computed. All topological combinations of duplicated (polyploid) nodes in the `multree` are compared, and the gene set with the lowest nRF value is used to generate `marker.list`.
- In `primer`, when selecting polyploid genomes, multi-FASTA files are generated for all gene copies of the target genes.