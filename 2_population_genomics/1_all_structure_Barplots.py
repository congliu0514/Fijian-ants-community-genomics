import ipyrad.analysis as ipa
import pandas as pd
import toyplot
import toyplot.pdf
import sys
import os
# 2. Barplots
for filename in os.listdir("/home/l/liu-cong/work/fijian_rad_project/4-analysis/structure/hdf_tem/"):
    #population map as dictioanry
    popfile=filename.split(".")[0]+"_popmap.txt"
    sp_name=filename.split(".")[0]
    imap = {}
    popmap=pd.read_table("/home/l/liu-cong/work/fijian_rad_project/4-analysis/structure/popmap/"+popfile,sep="\t")
    for i in popmap['pop'].unique().tolist():
        popmap_slice = popmap[popmap['pop'] == i]
        imap[i] = popmap_slice['indv'].tolist()
    
    struc = ipa.structure(
    data="/home/l/liu-cong/work/fijian_rad_project/4-analysis/structure/hdf_tem/"+filename,
    name=sp_name,
    workdir="analysis-structure/",
    imap=imap,
    load_only=True,
    )
    etable = struc.get_evanno_table([1,2, 3, 4, 5, 6])
    etable
    for k in range(2,8):
        table = struc.get_clumpp_table(k)
        # sort list by columns
        table.sort_values(by=list(range(k)), inplace=True)

        # or, sort by a list of names (here taken from imap)
        #import itertools
        #onames = list(itertools.chain(*imap.values()))
        #table = table.loc[onames]
        # build barplot
        canvas = toyplot.Canvas(width=500, height=250)
        axes = canvas.cartesian(bounds=("10%", "90%", "10%", "45%"))
        axes.bars(table)

        # add labels to x-axis
        ticklabels = [i for i in table.index.tolist()]
        axes.x.ticks.locator = toyplot.locator.Explicit(labels=ticklabels)
        axes.x.ticks.labels.angle = -60
        axes.x.ticks.show = True
        axes.x.ticks.labels.offset = 10
        axes.x.ticks.labels.style = {"font-size": "12px"}
        toyplot.pdf.render(canvas, sp_name+'_structure_k'+str(k)+'.pdf')