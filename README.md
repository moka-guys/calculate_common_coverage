# calculate_common_coverage
This script is designed to run in a docker container, utilising tools packaged in the container, namely sambamba (to run outside of docker the path ro sambamba mmust be changed).

## Inputs
The script has the following inputs:
- -r = minimum coverage level
- -m (boolean) = merge_overlapping_mate_reads
- -q = minimum base quality score
- -b = /path/to/bedfile.bed
- -o = /path/to/outputfile
- -t = /path/to/sample1_markdup.bam
- -u = /path/to/sample2_markdup.bam 
- -v = /path/to/sample3_markdup.bam (optional)

NB the BAM indexes must be present in the same folder as the BAM files

## How the script works

1. Takes a sambamba BED file and assesses for any overlapping regions within a single gene. A assertion error is raised if overlaps are found.

2. For each base sambamba (currently v0.7.1) is used to calculate the depth at a single base for all given BAM files using the following options:
  - -c minimum mean coverage for output (0)
  - -L genomic region in form chr:beg-end
  - -m detect overlaps of mate reads and handle them on per-base basis (if given as argument)
  - -q only count bases with a base quality score greater or equal to this value (if given)
  - -F filter reads using the filter ('mapping_quality >= 20')

3. The sambamba output is parsed, and the % of bases with sufficient coverage (where a base is covered above the minimum read depth in all samples) is reported for each gene.

## Output
The script creates a file (named using the path provided as an -o argument)

Output files match those created by the existing coverage tool, so they can be imported into MOKA using existing processes.

The output file contains one row per gene and three fields per row:
 1. entrezgene ID
 2. percentage of gene covered at min depth
 3. average coverage(not used in this analysis so a NULL value is returned).