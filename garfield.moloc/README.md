# GARFIELD output used to create the moloc bed file
1. GARFIELD was used to identify regions of the genome to test for colocalization of GWAS SNPs and molQTLs
2. GARFIELD identifies genomic regions in which molQTLs are enriched for GWAS SNPs. These regions are are used as input for the colocalization analysis, +/- 1 Mb from the lead SNP in the enriched region.
3. GARFIELD input must contain the following columns (chr, pos, snp, pval, anno_overlap)
	
	a. chr (chromosome)
	b. pos (chromosome position)
	c. snp (chromosome:position)
	d. pval (pvalue of the lead SNP; optional)
	e. anno_overlap (number of overlapping molQTL; max is 1; optional)

Note: columns 4 (pval) and 5 (anno_overlap) are optional but five columns are required. If the user does not wish to include this information, then they can fill the columns with 'NULL' values.