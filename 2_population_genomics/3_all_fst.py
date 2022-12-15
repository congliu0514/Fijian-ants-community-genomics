#Fst analysis
#import packages
import numpy as np
import h5py
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import bcolz
import pandas as pd
import allel
import os
count=1
for filename in os.listdir("../fst_vcf/"):
    sp_name=filename.split("_")[0]+"_"+filename.split("_")[1]
    print (str(count)+'/73')
    print (sp_name)
    count+=1
    
    #import vcf files
    vcf=allel.read_vcf("../fst_vcf/"+filename)
    df_sample=pd.read_csv('../fst_pop/'+sp_name+'_popmap.txt',sep='\t',index_col='index')


    #get the genotype from vcf file
    genotype=allel.GenotypeChunkedArray(vcf['calldata/GT'])
    pop_a=df_sample.population.unique()
    fst_result=pd.DataFrame(columns=pop_a,index=pop_a)
    for i in range(len(df_sample.population.unique())):
        for j in range(i+1,len(df_sample.population.unique())):
            vcf=allel.read_vcf("../fst_vcf/"+filename)
            #get the genotype from vcf file
            genotype=allel.GenotypeChunkedArray(vcf['calldata/GT'])
            df_sample=pd.read_csv('../fst_pop/'+sp_name+'_popmap.txt',sep='\t',index_col='index')
            #get two populations from the popfile
            pop1=df_sample.population.unique()[i]
            pop2=df_sample.population.unique()[j]
            n_sample_pop1=np.count_nonzero(df_sample.population==pop1)
            n_sample_pop2=np.count_nonzero(df_sample.population==pop2)
            #print (pop1,n_sample_pop1,pop2,n_sample_pop2)
            subpops={
            pop1:df_sample[df_sample.population==pop1].index,
            pop2:df_sample[df_sample.population==pop2].index,
            }
            acs = genotype.count_alleles_subpops(subpops)
            acu=allel.AlleleCountsArray(acs[pop1][:]+acs[pop2][:])
            flt=acu.is_segregating() & (acu.max_allele()==1)
            #print ('retaining',np.count_nonzero(flt),'SNPs')
            ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt,axis=0)[:,:2])
            ac2 = allel.AlleleCountsArray(acs[pop2].compress(flt,axis=0)[:,:2])
            genotype = genotype.compress(flt,axis=0)
            pop1_idx=subpops[pop1]
            pop2_idx=subpops[pop2]
            a,b,c = allel.stats.weir_cockerham_fst(genotype,subpops=[pop1_idx,pop2_idx],max_allele=1)
            snp_fst_wc = (a/(a+b+c))[:,0]
            #Calculate average FST
            fst_wc,se_wc,vb_wc,_ = allel.stats.blockwise_weir_cockerham_fst(genotype,subpops=[pop1_idx,pop2_idx],blen=10,max_allele=1)
            #print ('%.04f +/_ %.04f (weir & cockerham)' % (fst_wc,se_wc))
            fst_result.at[pop1,pop2]=fst_wc
            fst_result.at[pop2,pop1]=fst_wc

    fst_result.to_csv(sp_name+"_fst_matrix.csv")
        