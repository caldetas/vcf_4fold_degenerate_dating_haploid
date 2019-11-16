# vcf_4fold_degenerate_dating_haploid
A script for dating the branching points of different organisms, using the 4fold degenerate dating method, analysing the mutations in the protein-coding sequences.  


## Mutation Rate
The mutation rate used is the mutation Rate of *A.thaliana*.  
1.3*(10**-8) *Substitutions/year*
  
To adapt the parsing of gene-name search "parsing" in the script.  
Use option "-h" or "--help" for further information.  


## MANUAL

	
	vcf_4fold_degenerate_dating.py        
			
		--transcript <maker.fa>        
		--gff <maker.gff>        
		--vcf <vcf_file>        
			
		--cores <int> *specify nr of threads, default=cpu_count-1*        
		--mincov <minimum_coverages> *opt,comma separated integers in order of sample appearance in vcf*        
		--single *opt, consider only single variants(default!)*        
		--multi *opt, consider only multi variants*        
		--all *opt, consider single and multi variants*        
		--save_filtered_vcfs *opt,safe filtered lowcut vcf's for each sample*        
		--ffdg_pos_output *opt,safe ffdg positions on contigs*        
			
		*search for filtering (default: mindiff >=1 & ratio <=0.05)*        
		*pre-filter vcf with vcftools recommended, only snp's*        
		*for parsing adaptions search'parsing' in script*        
		*in multivariant mode each variant is seen as a own snp*        
		*change mutation rate in script (search 'substitution' in script)*        
		default is A.thaliana 1.3*(10**-8)(mut/(bp*year))        
			
	install packages with anaconda
	
