# vcf_4fold_degenerate_dating_haploid
a script to determine the branching points of different organisms  
using the 4fold degenerate dating method  
analysing only the protein-coding sequences  
  
  to adapt the parsing of gene-name search "parsing" in the script


## MANUAL

    vcf_4fold_degenerate_dating.py            
	            
	--mincov <minimum_coverages> *opt,comma separated integers in order of sample appearance in vcf*            
	--single *opt, consider only single variants(default=all)*            
	--multi *opt, consider only multi variants(default=all)*            
	--print_histo_df *opt,safe histogram df for R*            
	--save_filtered_vcfs *opt,safe filtered lowcut vcf's for each sample*            
	--ref <ref_sample_name> *opt,specify which sample is ref,filters all snps for ref*            
	--new_ref *define one of the samples as new outgroup*            
	--transcript <maker.fa>            
	--gff <maker.gff>            
	--vcf <vcf_file>            
	            
	*search for filtering (default: diff >=200 & ratio <=0.05)*            
	*pre-filter vcf with vcftools recommended*            
	*the filtering for the reference sample is inverted,            
	 filtering for ref mapping and valid snp in any other sample*            
	*new_ref specified takes all the valid snp's of the sample            
	 and inverts the var/ref values of all other samples*            
	*for parsing adaptions search'parsing' in script*            
	*in multivariant mode each variant is seen as a own spn*            
	            
install packages: https://scipy.org/install.html            
ubuntu: sudo apt-get install python-numpy python-pandas   
