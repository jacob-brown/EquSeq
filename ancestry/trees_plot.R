# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-15
# Last Modified: 2020-06-16

# Desc: 

source("../dependancies/treemix-1.13/src/plotting_funcs.R")



pdf("tree.pdf", 8, 8)
plot_tree("tree.out")
#plot_resid("tree.out", "poporder")
dev.off()
system("open -a Skim.app tree.pdf")




#prefix="tree.out"
#pdf("tree.pdf", 8, 8)
#par(mfrow=c(2,3))
#for(edge in 0:5){
#  plot_tree(cex=0.8,paste0(prefix,".",edge))
#  title(paste(edge,"edges"))
#}
#dev.off()
#system("open -a Skim.app tree.pdf")


