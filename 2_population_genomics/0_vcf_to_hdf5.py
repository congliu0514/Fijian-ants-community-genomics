# converting any VCF file to this HDF5 format

import ipyrad.analysis as ipa
import pandas as pd
import os
import toyplot

for filename in os.listdir("../vcf_filtered/"):
    converter = ipa.vcf_to_hdf5(
    name=filename,
    data="../vcf_filtered/"+filename,
    ld_block_size=20000,)
    converter.run()
