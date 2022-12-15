#Run STRUCTURE

import ipyrad.analysis as ipa
import pandas as pd
import toyplot
import sys
import os

for filename in os.listdir("/home/l/liu-cong/work/fijian_rad_project/4-analysis/structure/hdf_tem"):
    #population map as dictioanry
    popfile=filename.split(".")[0]+"_popmap.txt"

    imap = {}
    popmap=pd.read_table("/home/l/liu-cong/work/fijian_rad_project/4-analysis/structure/popmap/"+popfile,sep="\t")
    for i in popmap['pop'].unique().tolist():
        popmap_slice = popmap[popmap['pop'] == i]
        imap[i] = popmap_slice['indv'].tolist()
    # require that 50% of samples have data in each group
    minmap = {i: 0.1 for i in imap}
    # init analysis object with input data and (optional) parameter options
    struct = ipa.structure(
        name=filename.split(".")[0],
        data="/home/l/liu-cong/work/fijian_rad_project/4-analysis/jupyter/analysis-vcf2hdf5/"+filename,
        imap=imap,
        minmap=minmap,
        mincov=0.5,
    )


    # require that 50% of samples have data in each group
    struct.mainparams.burnin = 10000
    struct.mainparams.numreps = 200000
    struct.run(nreps=5, kpop=[1,2, 3, 4, 5, 6, 7], auto=True)
