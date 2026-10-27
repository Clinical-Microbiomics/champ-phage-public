# champ-phage-public

## License
The code in this repository is source-available for inspection; no license is granted.

## Purpose
This repository contains key scripts used in CHAMP phage profiling as of 7 September 2026. 
These are provided to allow inspection of code used to derive CHAMP phage relative abundance profiles.

## Scope & limitations
- This repository is provided as-is, without support.
- These scripts are provided for transparency and reproducibility purposes and are not intended to function as a complete, standalone pipeline.

## Contents
- Snakefile_phage_rules_only: Snakemake rules to run KMA mapping, rarefaction, and subsequent scripts.
- kmaprofiling.py: Performs filtering of KMA output.
- votu_summary.R: Estimates relative abundance from filtered read counts.
- votu_downsizeTarget.py: Determines rarefaction target.
- votu_summary_rarefy.R: Estimates rarefied relative abundance from rarefied filtered read counts.

## Further information
- The CHAMP profiling algorithms and benchmarks are described in Pita et al. (2024). 
  [Journal](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2024.1425489/full) / 
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/39483755/)
- CHAMP profiling is available exclusively from [Cmbio](https://cmbio.io). 

Copyright © 2026 Cmbio. All rights reserved.
