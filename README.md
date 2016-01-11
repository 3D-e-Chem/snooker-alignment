[![Build Status](https://travis-ci.org/3D-e-Chem/snooker-alignment.svg)](https://travis-ci.org/3D-e-Chem/snooker-alignment)
[![Codacy Badge](https://api.codacy.com/project/badge/grade/eebbe8d8006748e694e6a1668dcad33a)](https://www.codacy.com/app/sverhoeven/snooker-alignment)

Scripts and config file to make an sequence alignment from a bunch of GPCR protein sequences

See [snooker_align.vsd](snooker_align.vsd) for workflow.

The original alignment script was made for https://dx.doi.org/10.1186/1471-2105-12-332

# Requirements

* NCBI Blast
* hmmer
* Clustalo
* Perl packages:
  * Bio::SeqIO
  * Text::LevenshteinXS

# Gapped alignment based on gpcrdb human swissprot alignment

Steps to get an alignment

1. Create numbering schema
2. Download human swissprot alignment csv from gpcrdb website
3. Convert csv to with only positions of numbering schema
4. Create fasta from csv
5. Run blast with query seed alignment against swissprot/trembl/ensembl
5.1 Make sure all seed sequences have been found
6. Retrieve sequences of ids
7. Make sequences unique within same species
8. Remove species with less than 100 sequences
9. Run per tm alignment script
10. Remove sequences less than 9aa different within same species
11. Remove species with less than 100 sequences
12. Make tree of sequences
13. Generate entropy file based on tree

See [runs.md](runs.md) for commands to perform the steps.
