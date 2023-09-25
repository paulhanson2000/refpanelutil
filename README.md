## refpanelutil
This is a simple package of a few utility functions for working with LD reference panels.

* `liftOver()`: Use UCSC's [`liftOver`](https://genome.ucsc.edu/cgi-bin/hgLiftOver) program to convert genomic coordinates between builds. Includes conveniences such as automatically downloading the required chainfiles.
* `bcftoolsPipableCmd()`: An incomplete wrapper over [`bcftools`](https://samtools.github.io/bcftools/bcftools.html). Generates a single command made of a chain of `bcftools` commands, which can be piped. 
* `vcfs2Gds()`: Built ontop of `bcftoolsPipableCmd()` and [`SeqArray`](https://bioconductor.org/packages/release/bioc/html/SeqArray.html)'s `seqVCF2GDS()`, this function allows filtering and conversion of VCF files to be done together without any intermediate files being written.
* `ldCalc()`: Takes a `SeqArray` GDS file handle and calculates an LD matrix. It is possible to perform allele switching, and filter by region(s), variant IDs, and/or sample IDs. 

Note: the `liftOver()` function will not work on Windows, as `liftOver` is a Linux- & Mac-only program.
