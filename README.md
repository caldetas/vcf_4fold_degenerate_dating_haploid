# vcf_4fold_degenerate_dating_haploid
A script for dating the branching points of different organisms, using the 4fold degenerate dating method, analysing the mutations in the protein-coding sequences.  

## Mutation Rate
The mutation rate used is the mutation Rate of *A.thaliana*.  
1.3*(10**-8) *Substitutions/year*
  
To adapt the parsing of gene-name search "parsing" in the script.  
Use option "-h" or "--help" for further information.  
The "--ref" changes the refernce sample mutations (usually ~0). The option counts all the non mutated 4fdg sites, in which any other sample has a mutation, as a snp for the reference sample (gives total snp's organism).  
The "--new_ref" option is very usefull. Set a new reference sample to get better dating-resolution in the new-ref's branch.  

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
	*in multivariant mode each variant is seen as a own snp*      
	            
install packages: https://scipy.org/install.html            
ubuntu: sudo apt-get install python-numpy python-pandas   
