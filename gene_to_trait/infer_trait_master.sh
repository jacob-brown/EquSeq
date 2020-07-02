#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-29
# Last Modified: 2020-06-30
# Desc: 



scp jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/ephemeral/gene_to_trait/gl.out.* results/gene_to_trait/

gunzip  results/gene_to_trait/gl.out.beagle.gz

# convert beagle to non binary
python3 scripts/beagleBinary2Non.py results/gene_to_trait/gl.out.beagle results/gene_to_trait/gl.out.rn.beagle


# run the analysis
rm results/gene_to_trait/{beaglesub.csv,heat.pdf,heat.png,heat.comb.pdf,heat.comb.png}

python3 gene_to_trait/infer_trait.py 
Rscript gene_to_trait/simple_infer_trait.R
#Rscript gene_to_trait/infer_trait.R 