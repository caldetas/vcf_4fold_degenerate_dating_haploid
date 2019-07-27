# vcf_4fold_degenerate_dating_haploid
A script for dating the branching points of different organisms, using the 4fold degenerate dating method, analysing the mutations in the protein-coding sequences.  

## Mutation Rate
The mutation rate used is the mutation Rate of *A.thaliana*.  
1.3*(10**-8) *Substitutions/year*
  
To adapt the parsing of gene-name search "parsing" in the script.  
Use option "-h" or "--help" for further information.  
The "--ref" specifies the reference sample. The option helps to not get a row of zeros at the ref sample.  


## MANUAL

    vcf_4fold_degenerate_dating.py            
	--mincov <minimum_coverages> *opt,comma separated integers in order of sample appearance in vcf*            
	--single *opt, consider only single variants(default!)*            
	--multi *opt, consider only multi variants*            
	--all *opt, consider single and multi variants*            
	--print_histo_df *opt,safe histogram df for R*            
	--save_filtered_vcfs *opt,safe filtered lowcut vcf's for each sample*            
	--ffdg_pos_output *opt,safe ffdg positions on contigs*            
	--div_time_df *opt,safe divergence time df for each sample*            
	--ref <ref_sample_name> *opt,specify which sample is ref,filters all snps for ref*            
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
	*in multivariant mode each variant is seen as a own snp*            
	*change mutation rate in script, default is A.thaliana*            
	1.3*(10**-8)(mut/(bp*year))            

	            
install packages: https://scipy.org/install.html            
ubuntu: sudo apt-get install python-numpy python-pandas  
conda:	conda install -c conda-forge pandas
