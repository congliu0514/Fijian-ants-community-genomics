import ipyrad.analysis as ipa
import pandas as pd
import toyplot
import toyplot.pdf
import sys
import os
# 1. Choosing K
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
    etable = struc.get_evanno_table([1, 2, 3, 4, 5, 6])
    etable
    # get canvas object and set size
    canvas_k = toyplot.Canvas(width=400, height=300)

    # plot the mean log probability of the models in red
    axes = canvas_k.cartesian(ylabel="estLnProbMean")
    axes.plot(etable.estLnProbMean * -1, color="darkred", marker="o")
    axes.y.spine.style = {"stroke": "darkred"}

    # plot delta K with its own scale bar of left side and in blue
    axes = axes.share("x", ylabel="deltaK", ymax=etable.deltaK.max() + etable.deltaK.max() * .25)
    axes.plot(etable.deltaK, color="steelblue", marker="o");
    axes.y.spine.style = {"stroke": "steelblue"}

    # set x labels
    axes.x.ticks.locator = toyplot.locator.Explicit(range(len(etable.index)), etable.index)
    axes.x.label.text = "K (N ancestral populations)"
    toyplot.pdf.render(canvas_k, sp_name+'_k.pdf')    