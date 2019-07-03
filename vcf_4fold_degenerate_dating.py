#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 14:34:57 2018

@author: hannes
#!/usr/bin/env python3
"""

import sys
import getopt



def main(argv):

    #instructions
    form='vcf_4fold_degenerate_dating.py\
            \n\t\
            \n\t--mincov <minimum_coverages> *opt,comma separated integers in order of sample appearance in vcf*\
            \n\t--single *opt, consider only single variants(default=all)*\
            \n\t--multi *opt, consider only multi variants(default=all)*\
            \n\t--print_histo_df *opt,safe histogram df for R*\
            \n\t--save_filtered_vcfs *opt,safe filtered lowcut vcf\'s for each sample*\
            \n\t--ref <ref_sample_name> *opt,specify which sample is ref,filters all snps for ref*\
            \n\t--new_ref *define one of the samples as new outgroup*\
            \n\t--transcript <maker.fa>\
            \n\t--gff <maker.gff>\
            \n\t--vcf <vcf_file>\
            \n\t\
            \n\t*search for filtering (default: diff >=200 & ratio <=0.05)*\
            \n\t*pre-filter vcf with vcftools recommended*\
            \n\t*the filtering for the reference sample is inverted,\
            \n\t filtering for ref mapping and valid snp in any other sample*\
            \n\t*new_ref specified takes all the valid snp\'s of the sample\
            \n\t and inverts the var/ref values of all other samples*\
            \n\t*for parsing adaptions search\'parsing\' in script*\
            \n\t*in multivariant mode each variant is seen as a own snp*\
            \n\t\
            \ninstall packages: https://scipy.org/install.html\
            \nubuntu: sudo apt-get install python-numpy python-pandas\
            \n\t'
    
    import math    
    import pandas as pd
    import numpy as np
    import re

    #values
    vcf = ''
    variant = ''
    ref = 'NA' #for vcf's without ref-sample mapped: comment this line
    df_histo = 0
    save_filtered_vcfs = 0
    mincov = 0
    new_ref = ''


    try:
       opts, args = getopt.getopt(argv,"h",['mincov=', 'transcript=', 'gff=', 'vcf=', 'ref=', 'vcf=', 'new_ref=', 'print_histo_df', 'save_filtered_vcfs', 'single', 'multi'])
    except getopt.GetoptError:
       print ('{}'.format(form))
       sys.exit()
    for opt, arg in opts:
       if opt == '-h' or opt == '-help' or opt == '--help':
          print ('{}'.format(form))
          sys.exit()
       elif opt == '--transcript':
          transcript = arg
       elif opt == '--gff':
          gff = arg       
       elif opt == '--vcf':
          vcf = arg
       elif opt == '--single':
          variant = 'single'
       elif opt == '--multi':
          variant = 'multi'
       elif opt == '--mincov':
           mincov = arg.split(sep=',')
           for i in range(len(mincov)):
               mincov[i] = int(mincov[i])
       elif opt == '--print_histo_df':
           df_histo = 1  
       elif opt == '--save_filtered_vcfs':
           save_filtered_vcfs = 1  
       elif opt == '--ref':
           ref = arg
       elif opt == '--new_ref':
           new_ref = arg
           variant = 'single'


    mindiff = 1
    minfrac = 0.05


    #which column does the sample start?
    sample_col = 4















    def get_variant_histograms(vcf, new_ref):
        df = []
        ass=open('{}'.format(vcf))
        header = 0
        cnt=0
        for line in ass:
            cnt+=1
            if cnt % 1000 == 0:
                print('workin on vcf line:', cnt)
            line = line.strip('\n')
            line = line.split()
            if header == 0:
                 match= re.search("^(#CHROM)", line[0])
                 if match:
                    temp_list = ['contig']
                    temp_list.append(line[1])
                    #transition or transversion?
                    temp_list.append('substitution')
                    temp_list.append('multi_var')
                    #temp_list.extend(line[9:])
                    for i in range(len(line[9:])):
                        temp = line[9+i]
                        temp_list.extend(['{}_ref'.format(temp), '{}_var'.format(temp)])
                        temp = ''
                    #print(temp_list)
                    df.append(temp_list)
                    header = temp_list
                    temp_list = []
                    #debug
                    if new_ref != '':
                        for i in [i for i,x in enumerate(header) if x == '{}_ref'.format(new_ref)]:
                            new_ref_ref = i
                        for i in [i for i,x in enumerate(header) if x == '{}_var'.format(new_ref)]:
                            new_ref_var = i                    
            else:
                #first contig 
                if line[0] != '#':
                    match= re.search("^(\S+)", line[0])
                    if match:
                        line[0] = match.group(1)
                        temp_list = [line[0]]
                        temp_list.append(line[1])
                        #transition or transversion?
                        if  line[3] == 'A' and line[4] == 'G' or \
                            line[3] == 'G' and line[4] == 'A' or \
                            line[3] == 'C' and line[4] == 'T' or \
                            line[3] == 'T' and line[4] == 'C':
                                temp_list.append('ti')
                        else:
                                temp_list.append('tv')
                        temp_list.append('no') #multi_var
                        temp_list.extend(line[9:])
#                        print(temp_list)
                        #parsing AD-allele position
                        for i in (i for i,x in enumerate(line[8].split(':')) if x == 'AD'):
                            AD_pos = i
                        #select values reference/variance mappings
                        for i in range(len(temp_list[sample_col:])):
                            if temp_list[i+sample_col][0] != '.':
                                temp_list[i+sample_col] = temp_list[i+sample_col].split(sep=':')[AD_pos] #parsing vcf AD allele cov.
                            else:
                                temp_list[i+sample_col] = 'NA'
                        #new_ref sample specified??
                        if new_ref != '':
                            if len(temp_list[sample_col].split(sep=',')) >= 2:            
                                #print single variants
                                if len(temp_list[sample_col].split(sep=',')) == 2:
                                    
                                    #separate values for reference / variants on diff. columns
                                    temp_final = temp_list.copy()
                                    temp_final.extend(['']*len(temp_list[sample_col:]))
                                    for i in range(len(temp_list[sample_col:])):
#                                        print(temp_final[sample_col+i*2])
                                        if temp_list[i+sample_col] != 'NA':
                                            temp_final[sample_col+i*2] = temp_list[sample_col+i].split(sep=',')[0]
                                            temp_final[sample_col+i*2+1] = temp_list[sample_col+i].split(sep=',')[1]
                                        else:
                                            temp_final[sample_col+i*2] = np.NaN
                                            temp_final[sample_col+i*2+1] = np.NaN 
                                    #only select filtered values for good new_ref snps
                                    if float(temp_final[new_ref_ref]) - float(temp_final[new_ref_var]) >= float(mindiff):
                                        if float(temp_final[new_ref_var])/float(temp_final[new_ref_ref]) <= float(minfrac):
                                            df.append(temp_final)
#                                            print('original:')    #debug                                           
#                                            print(temp_final) #debug
#                                            print()    #debug 
                                    #invert all other samples if new_ref has a valid snp
                                    elif float(temp_final[new_ref_var]) - float(temp_final[new_ref_ref]) >= float(mindiff):
                                        if float(temp_final[new_ref_ref])/float(temp_final[new_ref_var]) <= float(minfrac):
#                                            print(temp_line)    #debug 
                                            for i in range(len(header[sample_col:])//2):
#                                                if header[sample_col +2*i][:-4] != new_ref:
                                                temp = temp_final[sample_col +2*i]
                                                temp_final[sample_col +2*i] = temp_final[sample_col +2*i +1]
                                                temp_final[sample_col +2*i +1] = temp
                                                temp = ''
                                            df.append(temp_final)
#                                            print('inverted:')    #debug                                           
#                                            print(temp_final) #debug
#                                            print()    #debug 
                        else:
                            if len(temp_list[sample_col].split(sep=',')) >= 2:            
                                #is it multi-var?
                                if len(temp_list[sample_col].split(sep=',')) > 2:
                                    for j in range(len(temp_list[sample_col].split(sep=','))-1):
                                        #create list
                                        temp_multi = temp_list.copy()
                                        #separate values for reference / variants on diff. columns
                                        temp_final = temp_multi.copy()
                                        temp_final.extend(['']*len(temp_multi[sample_col:]))
                                        #multi-var
                                        temp_final[3] = 'yes'
                                        #go trhough samples, make reference and variance column
                                        for i in range(len(temp_list[sample_col:])):
                                            if temp_list[i+sample_col] != 'NA':
                                                #separate values for reference / variants on diff. columns
                                                temp_final[sample_col+i*2] = temp_multi[sample_col+i].split(sep=',')[0]
                                                temp_final[sample_col+i*2+1] = temp_multi[sample_col+i].split(sep=',')[j+1]
                                                
                                            else:
                                                #write NA's in two columns
                                                temp_final[sample_col+i*2] = np.NaN #'NA'
                                                temp_final[sample_col+i*2+1] = np.NaN #'NA'                         
                                        #print line                                
                                        df.append(temp_final)
                                        temp_final = []
                                        temp_multi = []
                                #print single variants
                                elif len(temp_list[sample_col].split(sep=',')) == 2:
                                    #separate values for reference / variants on diff. columns
                                    temp_final = temp_list.copy()
                                    temp_final.extend(['']*len(temp_list[sample_col:]))
                                    for i in range(len(temp_list[sample_col:])):
#                                        print(temp_final[sample_col+i*2])
                                        if temp_list[i+sample_col] != 'NA':
                                            temp_final[sample_col+i*2] = temp_list[sample_col+i].split(sep=',')[0]
                                            temp_final[sample_col+i*2+1] = temp_list[sample_col+i].split(sep=',')[1]
                                        else:
                                            temp_final[sample_col+i*2] = np.NaN
                                            temp_final[sample_col+i*2+1] = np.NaN                      
                                    df.append(temp_final)
            temp_list = []

        ass.close()
        df = pd.DataFrame(df)
        df.columns = df.iloc[0]
        df = df.reindex(df.index.drop(0))
        df = df.reset_index(drop=True)
        for column in list(df.columns)[sample_col:]:
            df[column] = df[column].astype(np.float64)
        df['POS'] = df['POS'].astype(int)
        vcf_total = df
        del df
#        print('vcf_total') #debug
#        print(vcf_total) #debug
        return vcf_total
































    def filter_vcf(vcf, variant, ref, save_filtered_vcfs, df_histo, mincov):#, mindiff, minfrac):

        if ref == '':
            print()
            print('ERROR: ref-sample not specified')
            print()
            print ('{}'.format(form))
            sys.exit()        
        
        
        #filtering
#        global mindiff #diff from var coverage to ref coverage
#        global minfrac #fraction from variance coverage to reference coverage


        
        #make coverage tables with single/multi

#        print('get_variant_histograms(vcf)') #debug
        print()
        print()

        df = get_variant_histograms(vcf, new_ref)
        global df_vcf_total
        df_vcf_total = df

#        print(df) #debug


        lyst = list(df.columns)[sample_col:]
#        print(lyst) #debug
        #only single vcf

        if variant == 'single':
            df = df.loc[df.loc[:,'multi_var']== 'no']
        if variant == 'multi':
            df = df.loc[df.loc[:,'multi_var']== 'yes']
        
        if mincov != 0:
            if len(lyst)//2 != len(mincov):
                print()
                print()
                print('ERROR: number of mincov values does not match number of samples')
                print()
                print ('{}'.format(form))
                sys.exit()
    
        #df dictionary
        d_df = {}
    
#        print(df) #debug
        #variants calculation
        for i in range(len(lyst)//2):


#            print('type:', type(df.loc[0,lyst[2*i+1]]))
            #reference
            if lyst[2*i+1][:-4] == ref:
                    d_df[ref] = df.loc[   (df.loc[:,lyst[2*i]]-df.loc[:,lyst[2*i+1]] >= mindiff)   &   (df.loc[:,lyst[2*i+1]]/df.loc[:,lyst[2*i]] <= minfrac)   ]
                    if mincov != 0:
                        d_df[ref] = d_df[ref].loc[   d_df[ref].loc[:,lyst[2*i]] >= mincov[i]] 
            
            #variant
            else:
                    d_df[lyst[2*i+1][:-4]] = df.loc[   (df.loc[:,lyst[2*i+1]]-df.loc[:,lyst[2*i]] >= mindiff)   &   (df.loc[:,lyst[2*i]]/df.loc[:,lyst[2*i+1]] <= minfrac)   ]
                    if mincov != 0:
                        d_df[lyst[2*i+1][:]] = d_df[lyst[2*i+1][:]].loc[   d_df[lyst[2*i+1][:]].loc[:,lyst[2*i+1]] >= mincov[i]]
    
    
        #output on terminal (verbose)
        #safe df for return
        df_filtered_vcf = d_df.copy()
        #make smaller for ram
#        print(d_df)
        
    
        if mincov != 0:
            table = []
            for i in range(len(lyst)//2):
                table.append(['mincov {}:'.format(lyst[2*i+1][:])])
                table[-1].append(mincov[i])
            table = pd.DataFrame(table)
            
            table = table.set_index([0])
            table.columns = table.iloc[0]
            table = table.reindex(table.index.drop('mincov {}:'.format(lyst[2*0+1][:])))
            del table.index.name
            print()
            print()
            print('Minimal Coverages:')
            print()
            print(table)
    
    
    
        if df_histo == 1:
    
            pd.options.mode.chained_assignment = None  # default='warn' disabled
            df_rstudio = df.copy()
            df_rstudio['check'] = df_rstudio.contig.map(str) + ',' + df_rstudio.POS.map(str)
    
            #print(df)
            for i in list(d_df):
    
                d_df[i]['check'] = d_df[i].contig.map(str) + ',' + d_df[i].POS.map(str)
    
                
                df_rstudio.loc[~df_rstudio['check'].isin(d_df[i]['check']),('{}_ref'.format(i), '{}_var'.format(i))] = np.NAN
            del df_rstudio['check']
            df_rstudio.iloc[:,5:] = df_rstudio.iloc[:,5:].astype(float)
            df_rstudio = df_rstudio.fillna('NA')
            #save csv
            df_rstudio.to_csv(path_or_buf='df_rstudio_{}'.format(vcf\
                              .replace('.vcf', '')\
                              .replace('.gvcf', '')\
                              .replace('.g.vcf', '')), sep='\t', index = False)
    
            print(df_rstudio)
    
    
        if save_filtered_vcfs == 1:
    
    
            for i in range(len(lyst)//2):
                contig = ''
                vcf = open('{}'.format(vcf))
                out = open('{}.vcf'.format(lyst[i*2][:]), 'w')
                for line in vcf:
                    line = line.strip('\n')
                    #print commented lines
                    if line[0] == '#':
                        print(line,file=out)
                    else:
                        #set df of chr
                        if int(line.split()[0]) != contig:
                            contig = int(line.split()[0])
    
                            df_temp = d_df[lyst[i*2][:]].loc[   d_df[lyst[i*2][:]].loc[:,'contig'] == contig].loc[:,('contig','POS')]
                            df_temp = list(df_temp.loc[:,'POS'])
    
    
                        #print snp if in df
                        if int(line.split()[1]) in df_temp:
                            print(line,file=out)
    
                out.close()
                vcf.close()
    
        return df_filtered_vcf















    def get_pos_ffdg_on_contigs(transcript, gff):
        
        print()
        print('..reading transcript.fasta..')
        print()
        
        header = ''
        fa = open('{}'.format(transcript))
        seq = ''
        lyst_fa = {}
        first = 0
        for line in fa:
            line = line.strip('\n')
            if line[0] == '>':
                if first == 1:
                    lyst_fa[header] = seq
                    #print(lyst_fa[header])
                    header = line.split()[0][1:] #parsing fasta
                    #print(header)
                    
                    seq = ''
                else:
                    first = 1
                    header = line.split()[0][1:] #parsing fasta
                    #print(header)
    
            else:
                seq += line
                
        #last gene
        lyst_fa[header] = seq
        seq = ''
        fa.close()
#        print(lyst_fa) #debug
        print()
        print('..indexing fourfold-degenerate sites..')
        print()
        
    
        pos_fa = {}
        #4-fold degenerate codons
        fold = ['CT', 'GT', 'TC', 'CC', 'AC', 'GC', 'CG', 'GG']
        for gene in lyst_fa:
            pos_fa[gene] = []
#            print(pos_fa)
            for i in range(len(lyst_fa[gene])//3):
                if lyst_fa[gene][3*i:3*i+2] in fold:
                   pos_fa[gene].append(i*3+2+1)  #!!real position, +1
        del lyst_fa
        #print(pos_fa)
        
        print()
        print('..indexing gff gene positions from exons..')
        print()
        
        
        #open gff
        gff = open('{}'.format(gff))
        #make cds positions index
        cds_positions = {'check':[]}
        
        cc = {}
        for line in gff:
            if line[0] != '>':
                if line[0] != '#':
                    line = line.strip('\n')
                    line = line.split()
                    if len(line) >= 8:
                        if len(line[8].split(sep=';')) > 1:
    
                            if line[2] == 'exon':
#                                print('exon found') #debug
    
                                #parsing the gene name !!!important!!!
                                #has to be the same as the gene name parsed from cds fasta
                                #check if sep=';' or sep=':'
#                                gene = line[8].split(sep=';')[0].split(sep='=')[1].split('.')[0] #gatk
                                gene = line[8].split(sep=':')[0].split(sep='=')[1] #freebayes
                                if gene not in cc:
#                                    print(gene) #debug
                                    cc[gene] = []
                                    cc[gene].append([line[0], gene, line[3], line[4], line[6], 
                                       1+abs(int(line[3])-int(line[4]))]) #add 1! #parsing gff

                                else:
                                    cc[gene].append([line[0], gene, line[3], line[4], line[6], 
                                       1+abs(int(line[3])-int(line[4]))]) #add 1! #parsing gff
                                #parsing explained:
                                #[line[0], gene, line[3], line[4], line[6], 1+abs(int(line[3])-int(line[4]))]
                                #[contig, gene, start_bp, stop_bp, orientation(+/-), length(stop - start)]
                                
                                #append positions to cds_positions
                                for i in range(int(line[4])-int(line[3])+1):
                                    contig = line[0]
                                    pos = str(i + int(line[3]))
                                    cds_positions['check'].append(contig + ',' + pos)
            else:
                break
#                print('break')#debug
#        print('cc:', cc) #debug
        gff.close()
        cds_positions = pd.DataFrame(cds_positions)
#        print(cds_positions) #debug

        '''
        for i in cc:
            print(i)
            print(cc[i])       
        print('print(cc[i])') 
        '''                                                                                                                                     
        for i in cc:
            cc[i] = pd.DataFrame(cc[i])
            '''
            print(cc[i])
            '''
            #print(i)

            
#        print(cc)
#        print('cc')
#        print(pos_fa)
#        print('pos_fa')
        print()
        print('..calculating fourfold-degenerate positions on contigs..')
        print()
        '''
        for gene in pos_fa:
            if gene not in cc:
                print(gene)
        '''
        pos_con = [['contig', 'POS']]

        for gene in pos_fa:
#            print(gene)#debug
            #print(pos_fa[gene])
            #exception handling: exon has no snps
            if gene in cc:
#                print(cc[gene])
#                print('line164')
                for numb in pos_fa[gene]:
                    #print(pos_fa[gene])
                    numb = int(numb)
                    tot = 0
                    for i in range(len(cc[gene])):
                        tot += int(cc[gene].iloc[i,5])
#                        print(int(cc[gene].iloc[i,5]))
                        if tot >= numb:
                            diff = tot - numb
                            contig = str(cc[gene].iloc[i,0])
                            if cc[gene].iloc[i,4] == '+':
                                pos = int(cc[gene].iloc[i,3]) - diff
                            else:
                                pos = int(cc[gene].iloc[i,2]) + diff
                            pos_con.append([str(contig), int(pos)])
                            contig=''
                            pos = ''
                            break
#        print(pos_con) #debug
#        print('pos_con', len(pos_con)) #debug
        del pos_fa
        
        print()
        print('..determine 4fold-degenarate snps in snp.vcf..')
        print()

        checklist = pd.DataFrame(pos_con[1:], columns = pos_con[0])
#        print(checklist) #debug
#        sys.exit()
        del pos_con

        ffdg_positions_on_ctgs = checklist
#        print('ffdg_positions_on_ctgs', len(ffdg_positions_on_ctgs)) #debug
#        print(ffdg_positions_on_ctgs)
        return ffdg_positions_on_ctgs, cds_positions


















    def get_snp(transcript, gff, df_filtered_vcf, df_vcf_total, ffdg_positions_on_ctgs, cds_positions):
        filtered_vcf = df_filtered_vcf.astype(str) #necessary?
#        df = ffdg_positions_on_ctgs


#        print(df) #debug
 
        print()
        print('..compare snp.vcf with 4fold-degenarate candidates..')
        print()
       


        ffdg = len(ffdg_positions_on_ctgs) #total 4-fold degenerate sites

        if ffdg == 0:
            print()
            print('ERROR: no 4fold-dgenerate sites could be read from fasta')
            print()
            print(form)
            sys.exit()
        '''
        print(ffdg_positions_on_ctgs) #debug

        ffdg_positions_on_ctgs['check'] = ffdg_positions_on_ctgs.contig.map(str) + ',' + ffdg_positions_on_ctgs.POS.map(str)
        '''
        filtered_vcf['check'] = filtered_vcf.contig.map(str) + ',' + filtered_vcf.POS.map(str)
        df_ffdg_snps = filtered_vcf.loc[filtered_vcf['check'].isin(ffdg_positions_on_ctgs['check'])]
        df_cds_snps = filtered_vcf.loc[filtered_vcf['check'].isin(cds_positions['check'])]

#        print(df_ffdg_snps) #debug
#        print(filtered_vcf) #debug
#        print(cds_positions) #debug        
#        print(df_cds_snps) #debug 
        SNP = len(df_ffdg_snps)
        ti = len(df_ffdg_snps.loc[df_ffdg_snps.loc[:,'substitution'] == "ti"])
        tv = len(df_ffdg_snps.loc[df_ffdg_snps.loc[:,'substitution'] == "tv"])
        total = len(df_cds_snps)

#        print('SNP', SNP) #debug
#        print('ti', ti) #debug
#        print('tv', tv) #debug
#        print('ffdg', ffdg) #debug     
#        print('total', total) #debug   

        #calculate years apart with 4fdg-mutation constant
        #take in count the morigan-formula
        #weighting transitions and transversions

        #the K2P formula
        k2p = -0.5*math.log((1-2*(ti/ffdg)-(tv/ffdg)) * ((1-2*(tv/ffdg)))**.5)
        
        
        #the STD of K2p
        #split the whole moster-formula in four parts for easier editing
        t1 = 1/(1-2*(ti/ffdg)-(tv/ffdg))**2*(ti/ffdg)
        t2 = (0.5 * ((1/((1-2*(ti/ffdg)-(tv/ffdg))))+(1/(1-2*(tv/ffdg)))))**2*(tv/ffdg)
        t3 = 1/(1-2*(ti/ffdg)-(tv/ffdg))*(ti/ffdg);
        t4 = (0.5*((1/((1-2*(ti/ffdg)-(tv/ffdg))))+(1/(1-2*(tv/ffdg)))))*(tv/ffdg)
        
        k2p_std = ((1/ffdg)*(t1 + t2 -(t3 + t4)**2))**.5

        subst_rate = 1.3*(10**-8) #arabidopsis
#        subst_rate = 0.9*(10**-9) #asco 2002 min
#        subst_rate = 16.7*(10**-9) #asco 2002 max
        div_time     = k2p / (2*(subst_rate))
        
        div_time_std = k2p_std / (2*(subst_rate))

        if total != 0:
            ffdg_to_tot = SNP/total
        else:
            ffdg_to_tot = 'NA'
        #make stats df:
        stats = [   ["4fdg_sites:"   ,   ffdg],\
                    ["4fdg_snp:"   ,   SNP],\
#                    ["REF:"   ,   REF],\ 
                    ["snp_ratio:"   ,   SNP/ffdg],\
                    ["div_time"   ,   div_time],\
                    ["time_std"   ,   abs(div_time_std)],\
#                    ["transitions"   ,   ti],\
#                    ["transversions"   ,   tv],\
                    ["snp_tot:"   ,   total],\
                    ["4fdg/tot"   ,   ffdg_to_tot]\
                 ]

        print(pd.DataFrame(stats))
        print()
        print()
        #print('output: df_snp_{}'.format(vcf.replace('.vcf', '')))
        #print('output: df_snp_{}'.format(vcf.replace('.vcf', '')))
        #output.to_csv(path_or_buf='{}'.format('df_snp_{}'.format(vcf.replace('.vcf', ''))), sep='\t', index = False, header = False)
        #stats.to_csv(path_or_buf='{}'.format('stats_{}.txt'.format(vcf.replace('.vcf', ''))), sep='\t', index = False, header = False)
        

        return df_ffdg_snps, stats

































    df_filtered_vcf = filter_vcf(vcf, variant, ref, save_filtered_vcfs, df_histo, mincov)          #, mindiff, minfrac)
    ffdg_positions_on_ctgs, cds_positions = get_pos_ffdg_on_contigs(transcript, gff)
    '''
    print('print(df_filtered_vcf)')
    print(df_filtered_vcf)
    '''
    
#    for i in df_filtered_vcf:
#        print(i)
#    for i in df_filtered_vcf:
#        df_filtered_vcf[i] = df_filtered_vcf[i].iloc[:,:2]
#    print(df_filtered_vcf)

    print()
    print()
    print('***')
    print('vcf:',  vcf)
    print('***')
    print()
    print()

    #add 'check' column for comparing while get_snp()
#    print(ffdg_positions_on_ctgs) #debug
    ffdg_positions_on_ctgs['check'] = ffdg_positions_on_ctgs.iloc[:,0].map(str) + ',' + ffdg_positions_on_ctgs.iloc[:,1].map(str)


    results = {}
    #sample counter
    cnt = 0
    for i in df_filtered_vcf:
        cnt += 1
        print()
        print('calculating 4-fold snps of {}'.format(i), '\t', 'sample', cnt, 'of', len(df_filtered_vcf))
        print()
        output, stats = get_snp(transcript, gff, df_filtered_vcf[i], df_vcf_total, ffdg_positions_on_ctgs, cds_positions)
        results[i] = [output, stats]

  
    #make an index of samples
    lyst = []
    for i in df_filtered_vcf:
        lyst.append(i)

    #output formatting
    table = [['']]
    for j in range(len(results[lyst[0]][1])):
        table[0].append(results[lyst[0]][1][j][0])
    
    for i in lyst:
        table.append([i])
        for j in range(len(results[i][1])):
            table[-1].append(results[i][1][j][-1])


    df = pd.DataFrame(table)
    df.to_csv(path_or_buf='{}'.format('results.txt'), sep='\t', index = False, header = False)
    df = df.set_index([0])
    df.columns = df.iloc[0]
    df = df.reindex(df.index.drop(''))
    del df.index.name
    print()
    print()
    print()
    print('SUMMARY 4-fold degenerate sites:')
    print()
    print()
    #print whole results df
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(df)
    print()
    print()
    print()
    print()
    print('output:\tresults.txt')
    
    pd.options.mode.chained_assignment = None  # default='warn' disabled
    #same snps?
    '''
    for i in results:
        results[i][0] = results[i][0][results[i][0].iloc[:,2] == 'SNP']
        results[i][0]['check'] = results[i][0].iloc[:,0].map(str) + ',' + results[i][0].iloc[:,1].map(str)
    '''
    table = [['']]
    for x in lyst:
        table[0].append(x)
        table.append([x])
    
    for y in range(len(lyst)):
        for x in range(len(lyst)):
            table[y+1].append(str(len(results[lyst[x]][0].loc[results[lyst[x]][0]['check'] .isin( results[lyst[y]][0]['check'])])))

    df = pd.DataFrame(table)
    df.to_csv(path_or_buf='{}'.format('results_ffdg_snps_shared.txt'), sep='\t', index = False, header = False)
    df = df.set_index([0])
    df.columns = df.iloc[0]

    df = df.reindex(df.index.drop(''))
    del df.index.name
    print()
    print()
    print()
    print('SUMMARY shared 4-fold degenerate sites:')
    print()
    print()
    #print whole results df
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(df)
    print()
    print()
    print()
    print()
    print('output:\tresults_ffdg_snps_shared.txt')

    sys.exit()
    
if __name__ == "__main__":
    main(sys.argv[1:])
