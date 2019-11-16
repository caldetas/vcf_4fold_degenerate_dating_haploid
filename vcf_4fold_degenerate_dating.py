#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 14:34:57 2018

@author: hannes
#!/usr/bin/env python3
"""


import sys
import getopt
import math    
import pandas as pd
import numpy as np
import re
from multiprocessing import Pool,cpu_count


class data:
    
    def __init__(self):
            self.df = {}
            self.settings = {}

#init class data
data = data()


#instructions
form='\nvcf_4fold_degenerate_dating.py\
        \n\t\
        \n\t--transcript <maker.fa>\
        \n\t--gff <maker.gff>\
        \n\t--vcf <vcf_file>\
        \n\t\
        \n\t--cores <int> *specify nr of threads, default=cpu_count-1*\
        \n\t--mincov <int> *opt, one integer, multiple can be given for each sample comma separated in sample order of vcf columns*\
        \n\t--single *opt, consider only single variants(default!)*\
        \n\t--multi *opt, consider only multi variants*\
        \n\t--all *opt, consider single and multi variants*\
        \n\t--save_filtered_vcfs *opt,safe vcf with valid ffdg snp\'s for each sample pair*\
        \n\t--ffdg_pos_output *opt,safe ffdg positions on contigs*\
        \n\t\
        \n\t*search for filtering (default: mindiff >=1 & ratio <=0.05)*\
        \n\t*pre-filter vcf with vcftools recommended, only snp\'s*\
        \n\t*for parsing adaptions search\'parsing\' in script*\
        \n\t*in multivariant mode each variant is seen as a own snp*\
        \n\t*change mutation rate in script (search \'substitution\' in script)*\
        \n\tdefault is A.thaliana 1.3*(10**-8)(mut/(bp*year))\
        \n\t\
        \ninstall packages with anaconda\
        \n\t'


def prepare_multiprocess():
    
    """prepare data for multi-threading
    """
    
    data.settings['core_jobs'] = {}
    for i in range(data.settings['cores']):
        data.settings['core_jobs'][i] = []
    cnt=0
    for sample in data.settings['prefix']:
        data.settings['core_jobs'][(cnt)%data.settings['cores']].append(sample)
        cnt+=1
    data.settings['job_nr'] = list(range(len(data.settings['prefix'])))
    return


def pd_print(df):

    """Small function for printing pandas dataframe with header
    """

    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        df = df.reindex()
        df = df.set_index(df.iloc[:,0])
        df = df.reindex()
        del df[df.columns[0]]
        del df.index.name
        print(df)
    return


def filter_vcf(cnt):
    
    """function for reading the vcf file, setting a sample as reference
       and filtering for snps in all other samples
    """

    df_filtered_dic = {}
    sample_new_ref = data.settings['prefix'][cnt]
    cnt +=1
    df = data.df['vcf_original']#.copy()



        
    #reference searched for var
    df_var                      = df.loc[   (df.loc[:,'{}_var'.format(sample_new_ref)] - df.loc[:,'{}_ref'.format(sample_new_ref)]  >= data.settings['mindiff'])\
                                          & (df.loc[:,'{}_ref'.format(sample_new_ref)] / df.loc[:,'{}_var'.format(sample_new_ref)]  <= data.settings['minfrac'])\
                                                                        ]
    #reference searched for ref
    df_ref                      = df.loc[  (df.loc[:,'{}_ref'.format(sample_new_ref)] - df.loc[:,'{}_var'.format(sample_new_ref)]  >= data.settings['mindiff'])\
                                         & (df.loc[:,'{}_var'.format(sample_new_ref)] / df.loc[:,'{}_ref'.format(sample_new_ref)]  <= data.settings['minfrac'])\
                                                                        ]
    ref=[]
    if len(df_var) < 10:
        ref.append(sample_new_ref)
    if data.settings['mincov'] != 0:
        df_var = df_var.loc[df_var.loc[:,'{}_var'.format(sample_new_ref)] >= data.settings['mincov'][cnt-1]]
        df_ref = df_ref.loc[df.loc[:,'{}_ref'.format(sample_new_ref)] >= data.settings['mincov'][cnt-1]]
# =============================================================================
#             !!! ACHTUNG !!!
# =============================================================================
    
    for i, sample in enumerate(data.settings['prefix']):

                                                     #ref  compared with               vars
        df_filtered_vcf_ref = df_ref.loc [      (df_ref.loc[:,'{}_var'.format(sample)] - df_ref.loc[:,'{}_ref'.format(sample)]  >= data.settings['mindiff'])\
                                                                     & (df_ref.loc[:,'{}_ref'.format(sample)] / df_ref.loc[:,'{}_var'.format(sample)]  <= data.settings['minfrac'])\
                                                                        ]

                                                     #var  compared with               refs  
        df_filtered_vcf_var = df_var.loc [      (df_var.loc[:,'{}_ref'.format(sample)] - df_var.loc[:,'{}_var'.format(sample)]  >= data.settings['mindiff'])\
                                                                     & (df_var.loc[:,'{}_var'.format(sample)] / df_var.loc[:,'{}_ref'.format(sample)]  <= data.settings['minfrac'])\
                                                                        ]

        #more mincov filtering
        if data.settings['mincov'] != 0:
            df_filtered_vcf_ref = df_filtered_vcf_ref.loc[df_filtered_vcf_ref.loc[:,'{}_var'.format(sample)] >= data.settings['mincov'][i]]
            df_filtered_vcf_var = df_filtered_vcf_var.loc[df_filtered_vcf_var.loc[:,'{}_ref'.format(sample)] >= data.settings['mincov'][i]]

        #reduce to sliced df for storage
        df_filtered_vcf_ref = df_filtered_vcf_ref.loc[:,['contig', 'POS', 'substitution', 'check']]
        df_filtered_vcf_var = df_filtered_vcf_var.loc[:,['contig', 'POS', 'substitution', 'check']]

        #safe    
        df_filtered = pd.concat([df_filtered_vcf_ref, df_filtered_vcf_var])
        del df_filtered_vcf_ref
        del df_filtered_vcf_var
        df_filtered = df_filtered.sort_values(by=['contig', 'POS'], ascending=True)
        df_filtered_dic[sample] = df_filtered


        #print filtered vcf's
        if data.settings['save_filtered_vcfs'] == 1:

            contig = ''
            vcf = open('{}'.format(data.settings['vcf']))
            out = open('vcf_only_valid_ffdg_snps_{}_and_{}.vcf'.format(sample_new_ref, data.settings['prefix'][i]), 'w')
            for line in vcf:
                line = line.strip('\n')
                #print commented lines
                if line[0] == '#':
                    print(line,file=out)
                else:
                    #set df of chr
                    if int(line.split()[0]) != contig:
                        contig = int(line.split()[0])
                        df_temp = df_filtered_dic[data.settings['prefix'][i]].loc[   df_filtered_dic[data.settings['prefix'][i]].loc[:,'contig'] == str(contig)].loc[:,'POS'].values.tolist()
                   #print snp if in df
                    if int(line.split()[1]) in df_temp:
                        print(line,file=out)
            out.close()
            vcf.close()
    
    return df_filtered_dic,ref

def get_snp(cnt_sample):
    
    """function to to calculate divergence times with the substitution rate and the defined ffdg mutations
    """

    results = {}
    sample = data.settings['prefix'][cnt_sample]
    df_filtered_vcf = data.df['df_filtered_vcf'][sample]

    
    #sample counter
    cnt = 0

    for i in data.settings['prefix']:
        cnt += 1
        print('\n***\nworking on run {} of {} runs\n        on sample {} of {} samples\n***\n'.format(cnt, len(data.settings['prefix']), cnt_sample+1, len(data.settings['prefix'])))
        print()
        print('..compare snp.vcf with 4fold-degenarate candidates..')
        print()

        ffdg = len(data.df['ffdg_positions_on_ctgs']) #total 4-fold degenerate sites
        if ffdg == 0:
            print()
            print('ERROR: no 4fold-dgenerate sites could be read from fasta')
            print()
            print(form)
            sys.exit()

        df_ffdg_snps = df_filtered_vcf[i].loc[df_filtered_vcf[i].loc[:,'check'].isin(data.df['ffdg_positions_on_ctgs'].loc[:,'check'])]
        SNP = len(df_ffdg_snps)
        ti = len(df_ffdg_snps.loc[df_ffdg_snps.loc[:,'substitution'] == "ti"])
        tv = len(df_ffdg_snps.loc[df_ffdg_snps.loc[:,'substitution'] == "tv"])
        total = len(df_filtered_vcf[i])

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

# =============================================================================
#             here the substitution rate is chosen!!!
# =============================================================================

        subst_rate = 1.3*(10**-8) #arabidopsis substitution rate
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
                    ["snp_ratio:"   ,   SNP/ffdg],\
                    ["div_time"   ,   div_time],\
                    ["time_std"   ,   abs(div_time_std)],\
                    ["snp_tot:"   ,   total],\
                    ["4fdg/tot"   ,   ffdg_to_tot]\
                 ]

        print(pd.DataFrame(stats))
        print()
        print()

        results[i] = [df_ffdg_snps, stats]
    return results 


def printout(cnt):

    """function to create final distance matrix from calculated divergence times
    """


    print()
    print()
    print('***')
    print('vcf:',  data.settings['vcf'])
    print('***')
    print()
    print()

    sample = data.settings['prefix'][cnt]
    df_filtered_vcf = data.df['df_filtered_vcf'][sample]

    results = data.df['results'][sample]

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
    df = df.reindex()
    temp = df.iloc[0,0]
    df = df.set_index(df.iloc[:,0])
    df.columns = df.iloc[0,:]
    df = df.reindex(df.index.drop(temp))

    del df[temp]
    del df.index.name
    #store divergence dataframe in dictionary


    print(df)


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

    return df


def output():
    
    """function fro writing the csv of the calculated matrices
    """
    

    table = [['']]
    table[0].extend(data.settings['prefix'])
    for i_y, y in enumerate(data.settings['prefix']):
        table.append([y])

        for i_x, x in enumerate(data.settings['prefix']):
            table[-1].append(data.df['div_dic'][y].loc[x,'div_time'])
    table = pd.DataFrame(table[1:], columns = table[0])

    table = [['']]
    table[0].extend(data.settings['prefix'])
    for i_y, y in enumerate(data.settings['prefix']):
        table.append([y])

        for i_x, x in enumerate(data.settings['prefix']):
            table[-1].append(data.df['div_dic'][y].loc[x,'div_time'])
    table = pd.DataFrame(table[1:], columns = table[0])

    table_std = [['']]
    table_std[0].extend(data.settings['prefix'])
    for i_y, y in enumerate(data.settings['prefix']):
        table_std.append([y])

        for i_x, x in enumerate(data.settings['prefix']):
            table_std[-1].append(data.df['div_dic'][y].loc[x,'time_std'])
    table_std = pd.DataFrame(table_std[1:], columns = table_std[0])


    print('\nresults:\n')

    pd_print(table)
    table.to_csv('results_divergence_time.txt', sep = '\t', index = False)
    table_std.to_csv('results_divergence_time_std.txt', sep = '\t', index = False)


    print('\nOUTPUT:\nresults_divergence_time.txt\nresults_divergence_time_std.txt\n')
    if len(data.settings['ref']) > 0:
        print('Reference sample:', data.settings['ref'])


    return


def read_vcf_to_memory():
    df = []
    ass=open('{}'.format(data.settings['vcf']))
    header = 0
    cnt=0
    for line in ass:
        cnt+=1
        if cnt % 1000 == 0:
            print('reading vcf line:', cnt)
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
                #make sample index
                prefix = []
                for i in range(len(line[9:])):
                    temp = line[9+i]
                    prefix.append(temp)
                    temp_list.extend(['{}_ref'.format(temp), '{}_var'.format(temp)])
                    temp = ''
                df.append(temp_list)
                header = temp_list
                temp_list = []
           
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
                    #parsing AD-allele position
                    for i in (i for i,x in enumerate(line[8].split(':')) if x == 'AD'):
                        AD_pos = i
                    #select values reference/variance mappings
                    for i in range(len(temp_list[data.settings['sample_col']:])):
                        if temp_list[i+data.settings['sample_col']][0] != '.':
                            temp_list[i+data.settings['sample_col']] = temp_list[i+data.settings['sample_col']].split(sep=':')[AD_pos] #parsing vcf AD allele cov.
                        else:
                            temp_list[i+data.settings['sample_col']] = 'NA'

                    if len(temp_list[data.settings['sample_col']].split(sep=',')) >= 2:            
                        #is it multi-var?
                        if len(temp_list[data.settings['sample_col']].split(sep=',')) > 2:
                            for j in range(len(temp_list[data.settings['sample_col']].split(sep=','))-1):
                                #create list
                                temp_multi = temp_list.copy()
                                #separate values for reference / variants on diff. columns
                                temp_final = temp_multi.copy()
                                temp_final.extend(['']*len(temp_multi[data.settings['sample_col']:]))
                                #multi-var
                                temp_final[3] = 'yes'
                                #go trhough samples, make reference and variance column
                                for i in range(len(temp_list[data.settings['sample_col']:])):
                                    if temp_list[i+data.settings['sample_col']] != 'NA':
                                        #separate values for reference / variants on diff. columns
                                        temp_final[data.settings['sample_col']+i*2] = temp_multi[data.settings['sample_col']+i].split(sep=',')[0]
                                        temp_final[data.settings['sample_col']+i*2+1] = temp_multi[data.settings['sample_col']+i].split(sep=',')[j+1]
                                        
                                    else:
                                        #write NA's in two columns
                                        temp_final[data.settings['sample_col']+i*2] = np.NaN #'NA'
                                        temp_final[data.settings['sample_col']+i*2+1] = np.NaN #'NA'                         
                                #print line                                
                                df.append(temp_final)
                                temp_final = []
                                temp_multi = []
                        #print single variants
                        elif len(temp_list[data.settings['sample_col']].split(sep=',')) == 2:
                            #separate values for reference / variants on diff. columns
                            temp_final = temp_list.copy()
                            temp_final.extend(['']*len(temp_list[data.settings['sample_col']:]))
                            for i in range(len(temp_list[data.settings['sample_col']:])):
                                if temp_list[i+data.settings['sample_col']] != 'NA':
                                    temp_final[data.settings['sample_col']+i*2] = temp_list[data.settings['sample_col']+i].split(sep=',')[0]
                                    temp_final[data.settings['sample_col']+i*2+1] = temp_list[data.settings['sample_col']+i].split(sep=',')[1]
                                else:
                                    temp_final[data.settings['sample_col']+i*2] = np.NaN
                                    temp_final[data.settings['sample_col']+i*2+1] = np.NaN                      
                            df.append(temp_final)
        temp_list = []
    #store the vcf table
    ass.close()
    df = pd.DataFrame(df)
    df.columns = df.iloc[0]
    df = df.reindex(df.index.drop(0))
    df = df.reset_index(drop=True)
    for column in list(df.columns)[data.settings['sample_col']:]:
        df[column] = df[column].astype(np.float64)
    df['POS'] = df['POS'].astype(int)
    df= df.copy()
    data.df['vcf_original'] = df
    data.df['vcf_original'].loc[:,'check'] = data.df['vcf_original'].contig.map(str) + ',' + data.df['vcf_original'].POS.map(str)
    #more filtering from commandline settings
    if data.settings['variant'] == 'single':
        data.df['vcf_original'] = data.df['vcf_original'].loc[data.df['vcf_original'].loc[:,'multi_var']== 'no']
    if data.settings['variant'] == 'multi':
        data.df['vcf_original'] = data.df['vcf_original'].loc[data.df['vcf_original'].loc[:,'multi_var']== 'yes']
    #store sample index
    data.settings['prefix'] = prefix
    if data.settings['mincov'] != 0:
        if len(data.settings['mincov']) == 1:
            data.settings['mincov'] = len(data.settings['prefix'])*data.settings['mincov']
        elif len(data.settings['mincov']) != len(data.settings['prefix']):
            print()
            print('ERROR: number of mincov values not equal with number of samples')
            print()
            sys.exit()
    return





def get_pos_ffdg_on_contigs():
    
    print()
    print()
    print('..reading transcript.fasta..')
    print()
    genes_contain_points = 'no'
    header = ''
    fa = open('{}'.format(data.settings['transcript']))
    seq = ''
    lyst_fa = {}
    first = 0
    for line in fa:
        line = line.strip('\n')
        if line[0] == '>':
            if first == 1:
                lyst_fa[header] = seq
                header = line.split()[0][1:] #parsing fasta header for gene name
                seq = ''
                if '.' in header:
                    genes_contain_points = 'yes'
            else:
                first = 1
                if '.' in header:
                    genes_contain_points = 'yes'
                header = line.split()[0][1:] #parsing fasta header for gene name
        else:
            seq += line
            
    #last gene
    lyst_fa[header] = seq
    seq = ''
    fa.close()
    print()
    print('..indexing fourfold-degenerate sites..')
    print()

    pos_fa = {}
    fold = ['CT', 'GT', 'TC', 'CC', 'AC', 'GC', 'CG', 'GG']
    for gene in lyst_fa:
        pos_fa[gene] = []
        for i in range(len(lyst_fa[gene])//3):
            if lyst_fa[gene][3*i:3*i+2] in fold:
               pos_fa[gene].append(i*3+2+1)  #!!real position, +1
    del lyst_fa

    print()
    print('..indexing gff gene positions from exons..')
    print()

    #open gff
    gff = open('{}'.format(data.settings['gff']))
    #make cds positions index
    cds_positions = [['contig', 'POS']]
    
    cc = {}
    for line in gff:
        if line[0] != '>':
            if line[0] != '#':
                line = line.strip('\n')
                line = line.split()
                if len(line) >= 8:
                    if len(line[8].split(sep=';')) > 1:
                        if line[2] == 'exon':
# =============================================================================
#                                 parsing the gene name !!!important!!!
#                                 has to be the same as the gene name parsed from cds fasta
#                                 check if sep=';' or sep=':'
# =============================================================================

                            #gene names contain points??
                            if genes_contain_points == 'no':
                                gene = line[8].split(sep='ID=')[1].split(sep=':')[0].split(sep=';')[0].split('.')[0] #Bgt
                            elif genes_contain_points == 'yes':
                               gene = line[8].split(sep='ID=')[1].split(sep=':')[0].split(sep=';')[0]                #Cladonia
                            if gene not in cc:
                                cc[gene] = []
                                cc[gene].append([line[0], gene, line[3], line[4], line[6], 
                                   1+abs(int(line[3])-int(line[4]))]) #add 1! #parsing gff

                            else:
                                cc[gene].append([line[0], gene, line[3], line[4], line[6], 
                                   1+abs(int(line[3])-int(line[4]))]) #add 1! #parsing gff
# =============================================================================
#                                 parsing explained:
#                                 [line[0], gene, line[3], line[4], line[6], 1+abs(int(line[3])-int(line[4]))]
#                                 [contig, gene, start_bp, stop_bp, orientation(+/-), length(stop - start)]
#                                 
#                                 append positions to cds_positions
# =============================================================================
                            for i in range(int(line[4])-int(line[3])+1):
                                contig = line[0]
                                pos = str(i + int(line[3]))
                                cds_positions.append([contig, pos])
        else:
            break
    gff.close()
    cds_positions = pd.DataFrame(cds_positions[1:], columns=cds_positions[0])
                                                                                                           
    for i in cc:
        cc[i] = pd.DataFrame(cc[i])


    print()
    print('..calculating fourfold-degenerate positions on contigs..')
    print()

    pos_con = [['contig', 'POS']]

    for gene in pos_fa:

        if gene in cc:

            for numb in pos_fa[gene]:

                numb = int(numb)
                tot = 0

                for i in range(len(cc[gene])):
                    tot += int(cc[gene].iloc[i,5])

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

    del pos_fa
    
    print()
    print('..determine 4fold-degenarate snps in snp.vcf..')
    print()



    ffdg_positions_on_ctgs = pd.DataFrame(pos_con[1:], columns = pos_con[0])
    del pos_con
    data.df['ffdg_positions_on_ctgs'] = ffdg_positions_on_ctgs
    data.df['cds_positions'] = cds_positions
    
    #store ffdgs positon in csv
    if  data.settings['ffdg_pos_output'] == 1:
         data.df['ffdg_positions_on_ctgs'].iloc[:,:2].to_csv('ffdg_positions_on_ctgs', sep='\t', index=False)
    data.df['ffdg_positions_on_ctgs'].loc[:,'check'] = data.df['ffdg_positions_on_ctgs'].iloc[:,0].map(str) + ',' + data.df['ffdg_positions_on_ctgs'].iloc[:,1].map(str)
#    data.df['cds_positions'].loc[:,'check'] = data.df['cds_positions'].iloc[:,0].map(str) + ',' + data.df['cds_positions'].iloc[:,1].map(str)
    return






if __name__ == "__main__": 

    #define variables in class
    data.settings['vcf'] = 0
    data.settings['variant'] = 'single'
    data.settings['save_filtered_vcfs'] = 0
    data.settings['mincov'] = 0
    data.settings['ffdg_pos_output'] = 0
    data.settings['cores'] = int(cpu_count()-1)
    data.settings['ref'] = [] #list with possible reference candidates
    data.df['results'] = {}
    data.df['df_filtered_vcf'] = {}
    data.df['div_dic'] = {}

    #filtering options
    data.settings['mindiff'] = 1
    data.settings['minfrac'] = 0.05


    #which column does the sample start in temp cov df?
    data.settings['sample_col'] = 4

    try:
       opts, args = getopt.getopt(sys.argv[1:],"h",['cores=', 'mincov=', 'transcript=', 'gff=', 'vcf=', 'ffdg_pos_output',  'save_filtered_vcfs', 'single', 'multi', 'all'])
    except getopt.GetoptError:
       print ('{}'.format(form))
       sys.exit()
    for opt, arg in opts:
       if opt == '-h' or opt == '-help' or opt == '--help':
          print ('{}'.format(form))
          sys.exit()
       elif opt == '--transcript':
          data.settings['transcript'] = arg
       elif opt == '--gff':
          data.settings['gff'] = arg
       elif opt == '--vcf':
           data.settings['vcf'] = arg
       elif opt == '--single':
          data.settings['variant'] = 'single'
       elif opt == '--all':
          data.settings['variant'] = ''
       elif opt == '--multi':
          data.settings['variant'] = 'multi'
       elif opt == '--mincov':
           mincov = arg.split(sep=',')
           for i in range(len(mincov)):
               mincov[i] = int(mincov[i])
               data.settings['mincov'] = mincov
       elif opt == '--print_histo_df':
           data.settings['print_histo_df'] = 1  
       elif opt == '--save_filtered_vcfs':
           data.settings['save_filtered_vcfs'] = 1 
       elif opt == '--ffdg_pos_output':
           data.settings['ffdg_pos_output'] = 1
       elif opt == '--div_time_df':
           data.settings['div_t_df'] = 1
       elif opt == '--cores':
           data.settings['cores'] = int(arg)
    if sys.argv[1:] == []:
          print ('{}'.format(form))
          sys.exit()







# =============================================================================
#     Execution of programs start here!!!
# =============================================================================


    read_vcf_to_memory()
    get_pos_ffdg_on_contigs()


    #prepare list with jobs
    prepare_multiprocess()


    #execute multi-threading
    p = Pool(data.settings['cores']) #number of processes
    filtered_vcf = p.map(filter_vcf, data.settings['job_nr']) 

    #store df to class
    for x,sample in enumerate(data.settings['prefix']):
        data.df['df_filtered_vcf'][sample] = filtered_vcf[x][0]
        if len(filtered_vcf[x][1]) > 0:
            data.settings['ref'].extend(filtered_vcf[x][1])
    data.settings['ref'] = set( data.settings['ref'])


    #execute multi-threading
    p = Pool(data.settings['cores'])
    snp = p.map(get_snp, data.settings['job_nr']) #execute multi-threading

    #store df to class
    for x,sample in enumerate(data.settings['prefix']):
        data.df['results'][sample] = snp[x]


    #execute multi-threading
    p = Pool(data.settings['cores'])
    printed = p.map(printout, data.settings['job_nr']) #execute multi-threading

    #store df to class
    for x,sample in enumerate(data.settings['prefix']):
        data.df['div_dic'][sample] = printed[x]


    #final output
    output()


    sys.exit()