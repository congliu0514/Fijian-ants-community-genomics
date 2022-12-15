#This script converts ipyrad vcf file to vcf with single snp per locus.
#usage: python ipyrad_vcf2single_snp.py vcf1 vcf2
import sys
print 'converts ipyrad vcf file vcf with single snp per locus'
vcf=open(sys.argv[1])
vcf2=open(sys.argv[2],'w')
n_sample=int(sys.argv[3])
seen = set() # set for fast O(1) amortized lookup
for line in vcf:
    if line.startswith('#'):
        vcf2.write(line)
    else:
        z =line.split('\t')
        if z[0] not in seen:
            seen.add(z[0])
            for i in range(n_sample+9):
                if z[i][-1]!='\n':
                    vcf2.write(z[i]+'\t')
                else:
                    vcf2.write(z[i])
vcf2.close()
