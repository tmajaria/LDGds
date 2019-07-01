# LDGds
## Calculate linkage disequilibrium from a GDS file.
Maintainer: Analysis Commons  
Version: 0.1

## Description:
Generate LD measures from genotypes in [Genomic Data Structure (GDS)](https://www.biostat.washington.edu/sites/default/files/modules/GDS_intro.pdf) format. This workflow will return LD information for a set of defined samples over a set of variants or a defined variant range. A flat file of LD values and a simple visualization are returned

### What data are required for this workflow to run?
This workflow requires genotype files in GDS format (\*.gds) and a reference variant or genomic region.

### What does this workflow output?
An M x N diagonal matrix of pairwise LD values where N = number of variants within the genomic region. M = 1 if a reference variant is provided. M = N otherwise.

## Workflow Inputs
**Bold** inputs are required. *Italic* inputs are optional. **NOTE**: one of [reference variant, interval] is required. Both may be provided resulting in LD values for the reference variant within the interval.

- **this_gds_file**: [file, \*.gds] GDS file of genotypes per sample.
- *this_sample_ids_file*: [file, default = all samples] File of sample IDs desired for LD calculation. This file should contain one sample ID per line with no header.
- *this_ref_var*: [string, chr:pos] Genetic variant for which LD should be calculated. If provided, output is a row vector with pairwise LD with this variant in each row entry. Variant format should be 'chromosome:position'. Any punctuation seperator may be used. Only the first two values separated by punctuation will be considered.
- *this_interval*: [string, chr:start:end] Genomic interval for whcih LD should be calculated. If provided, LD will be calculated for only those variants falling within this interval. Interval format should be 'chromosome:start:end'. Any punctuation seperators may be used and need not match. Only the first three values separated by punctuation will be considered.
- *this_half_interval*: [int, default = 25kb] 1/2 of desired interval length if no interval is provided. When only a reference variant is provided, this value will be added and subtracted from the reference variant position to define the interval end and start, respectively. 
- *this_min_mac*: [int, default = 0] Minimum minor allele count for variant to be included in LD calculation.
- *this_max_mac*: [int, default = inf] Maximum minor allele count for variant to be included in LD calculation.
- *this_min_maf*: [int, default = 5%] Minimum minor allele frequency for variant to be included in LD calculation.
- *this_max_maf*: [int, default = 1] Maximum minor allele frequency for variant to be included in LD calculation.
- *this_ld_method*: [string, default = 'r'] LD calculation method. This value refers to the output LD values. Refer to documentation for the [SNPRelate package](https://bioconductor.org/packages/release/bioc/html/SNPRelate.html), specifically the function *snpgdsLDMat*. Possible values are: composite, correlation, r, and dprime with reasonable abbreviations accepted.
- **this_out_pref**: [string] Prefix for output file.
- *this_memory*: [int, default = 5GB] Amount of memory to request for computation in GB.
- *this_disk*: [int, default = size(this_gds_file) + 20 GB] Amount of disk space to request for computation in GB.
  
## Workflow Output

- **ld_file**: [file, \*.csv] An M x N diagonal matrix of pairwise LD values where N = number of variants within the genomic region. M = 1 if a reference variant is provided. M = N otherwise.