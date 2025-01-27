# Flying-plankton
This project analyzes dispersal and colonization abilities in microorganisms from a trait perspective in order to understand the contribution of phenotypic and metabolic traits shaping microbial communities

# Flying-plankton
Shurin Lab (UCSD) flying plankton collaboration

#### Associated article:
Traits Determine Dispersal and Colonization Abilities of Microbes
#### Authors
Isidora Echenique-Subiabre, Sara L. Jackrel, Jay McCarren, Chase C. James, Elisabet Perez-Coronel, Cindy Tran, Madeline Perreault, Ugbad Farah, P. Signe White, Henry K. Baker, Christopher Wall, Lindsay Sager, Scott Becker, Andrew Barton and Jonathan B. Shurin

[![DOI](https://zenodo.org/badge/631060996.svg)](https://doi.org/10.5281/zenodo.14744466)


# Notes on cleaned data file

Hi all,

A few comments on the cleaned data I just uploaded.

If you haven't used .Rdata files before, use the function load() to pull in these files

This will load in 3 data frames directly to your environment: asv_table, taxonomy, and metadata

Quick description of each:

asv_table: data frame where rows are samples and columns are ASVs

taxonomy: data frame with the ASV Hash, full taxonomy, and taxonomy split by level (columns A-G)

metadata: data frame with metadata for samples, row names for asv_table align with the column "sample_name" within this data frame.
          rows should be in the same order as rows for the asv_table but please use the match() function or a similar function when linking the two data frames together.

