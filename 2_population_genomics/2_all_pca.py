# all pca

import ipyrad.analysis as ipa
import pandas as pd
import os
import toyplot
import toyplot.pdf

for filename in os.listdir("/home/l/liu-cong/work/fijian_rad_project/4-analysis/jupyter/analysis-vcf2hdf5/"):
    sp_name=filename.split(".")[0]
    #population map as dictioanry
    popfile=sp_name+"_popmap.txt"

    imap = {}
    popmap=pd.read_table("/home/l/liu-cong/work/fijian_rad_project/4-analysis/structure/popmap/"+popfile,sep="\t")
    for i in popmap['pop'].unique().tolist():
        popmap_slice = popmap[popmap['pop'] == i]
        imap[i] = popmap_slice['indv'].tolist()
    # require that 10% of samples have data in each group
    minmap = {i: 0.1 for i in imap}
    pca = ipa.pca(
    data=data,
    imap=imap,
    minmap=minmap,
    mincov=0.75,
    impute_method="sample",
    )
    # run the PCA analysis
    pca.run()
    # store the PC axes as a dataframe
    df = pd.DataFrame(pca.pcaxes[0], index=pca.names)
    # write the PC axes to a CSV file
    df.to_csv(sp_name+"_pca_analysis.csv")
    # save returned values of the plot command
    canvas, axes, mark = pca.run_and_plot_2D(0, 1, subsample=False)
    # pass the canvas object to render function
    toyplot.pdf.render(canvas, sp_name+"-PCA.pdf")
    pca.run_tsne(subsample=True, perplexity=4.0, n_iter=100000, seed=123)
    canvas, axes, mark =pca.draw();
    toyplot.pdf.render(canvas, sp_name+"-PCA_TNSE.pdf")