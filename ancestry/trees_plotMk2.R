# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-07-07
# Last Modified: 2020-07-07

# Desc: 

# cd sandbox/plotting_trees
# cp ../../results/ancestry/treemix/tree.out.1* .
###########################################
################# Modules #################
###########################################

library(dplyr)
library(ggplot2)

################################
### plot trees


source("../../scripts/jb_plotting_funcs.R") # adds a xmax = val

### single migrations ###
pdf("tree.pdf", 20, 15)
plot_tree("tree.out.1", plus=0.0001, disp=0.00003)
options(scipen = 0)
dev.off()
system("open -a Skim.app tree.pdf")



#plot_tree("tree.out.0", o = NA, cex = 1, disp = 0.003, plus = 0.01, flip = vector(), arrow = 0.05, scale = T, ybar = 0.1, mbar = T, plotmig = T, plotnames = T, xmin = 0, xmax = 0.0006,  lwd = 1, font = 1)

######
#
stem="tree.out.1"
o = NA
cex = 1
disp = 0.003
plus = 0.01
flip = vector()
arrow = 0.05
scale = T
ybar = 0.1
mbar = T
plotmig = T
plotnames = T
xmin = 0
xmax = F
lwd = 1
font = 1

d = paste(stem, ".vertices.gz", sep = "")
e = paste(stem, ".edges.gz", sep = "")
se = paste(stem, ".covse.gz", sep = "")
d = read.table(gzfile(d), as.is = T, comment.char = "", quote = "")
e = read.table(gzfile(e), as.is  = T, comment.char = "", quote = "")
if (!is.na(o)){
	o = read.table(o, as.is = T, comment.char = "", quote = "")
}
e[,3] = e[,3]*e[,4]
e[,3] = e[,3]*e[,4]

se = read.table(gzfile(se), as.is = T, comment.char = "", quote = "")
m1 = apply(se, 1, mean)
m = mean(m1)
#m = 0
for(i in 1:length(flip)){
	d = flip_node(d, flip[i])
}
d$x = "NA"
d$y = "NA"
d$ymin = "NA"
d$ymax = "NA"
d$x = as.numeric(d$x)
d$y = as.numeric(d$y)
d$ymin = as.numeric(d$ymin)
d$ymax = as.numeric(d$ymax)

d = set_y_coords(d)
d = set_x_coords(d, e)
print(d)
d = set_mig_coords(d, e)

pdf("tree.pdf", 20, 15)
plot(d$x, d$y, axes = F, ylab = "", xlab = "Drift parameter", xlim = c(xmin, max(d$x)+plus), pch = "")
axis(1)
mw = max(e[e[,5]=="MIG",4])
mcols = rev(heat.colors(150))
for(i in 1:nrow(e)){
	col = "black"
	if (e[i,5] == "MIG"){
		w = floor(e[i,4]*200)+50
		if (mw > 0.5){
			w = floor(e[i,4]*100)+50
		}
		col = mcols[w]
		if (is.na(col)){
			col = "blue"
		}
	}
	v1 = d[d[,1] == e[i,1],]
	v2 = d[d[,1] == e[i,2],]
	if (e[i,5] == "MIG"){
		if (plotmig){
		arrows( v1[1,]$x, v1[1,]$y, v2[1,]$x, v2[1,]$y, col = col, length = arrow)
		}
	}
	else{
		lines( c(v1[1,]$x, v2[1,]$x), c(v1[1,]$y, v2[1,]$y), col = col, lwd = lwd)  #### its something about the x-y lim values
	}
}


tmp = d[d[,5] == "TIP",]
print(tmp$x)
print(disp)
if ( !is.na(o)){
	for(i in 1:nrow(tmp)){
		tcol = o[o[,1] == tmp[i,2],2]
		if(plotnames){
			#print(tmp[i,2])
			text(tmp[i,]$x+disp, tmp[i,]$y, labels = tmp[i,2], adj = 0, cex = cex, col  = tcol, font = font)
		}
	}
}else{
	if (plotnames){
	text(tmp$x+disp, tmp$y, labels = tmp[,2], adj = 0, cex = cex, font = font)
	}
}
dev.off()
system("open -a Skim.app tree.pdf")
if (scale){
print (paste("mse", mse))
    lines(c(0, mse*10), c(ybar, ybar))
text( 0, ybar - 0.04, lab = "10 s.e.", adj = 0, cex  = 0.8)
lines( c(0, 0), c( ybar - 0.01, ybar+0.01))
lines( c(mse*10, mse*10), c(ybar- 0.01, ybar+ 0.01))
}
    if (mbar){
            mcols = rev( heat.colors(150) )
            mcols = mcols[50:length(mcols)]
            ymi = ybar+0.15
            yma = ybar+0.35
            l = 0.2
            w = l/100
            xma = max(d$x/20)
            rect( rep(0, 100), ymi+(0:99)*w, rep(xma, 100), ymi+(1:100)*w, col = mcols, border = mcols)
            text(xma+disp, ymi, lab = "0", adj = 0, cex = 0.7)
	if ( mw >0.5){ text(xma+disp, yma, lab = "1", adj = 0, cex = 0.7)}
            else{
		text(xma+disp, yma, lab = "0.5", adj = 0, cex =0.7)
	}
	text(0, yma+0.06, lab = "Migration", adj = 0 , cex = 0.6)
	text(0, yma+0.03, lab = "weight", adj = 0 , cex = 0.6)
    }	
}


###########################################
############## Function(s) ################
###########################################
require(ape)

tree_all <- read.tree(gzfile("tree.out.1.treeout.gz"))
tree <- tree_all[[1]]
str(tree)


#tree$edge <- read.table(gzfile("tree.out.0.edges.gz"), as.is  = T, comment.char = "", quote = "")

pdf("ape.pdf", 15, 15)
plot(tree)
dev.off()
system("open -a Skim.app ape.pdf")


pdf("ape.pdf", 15, 15)
plot.phylo(tree,use.edge.length = T,no.margin = T,type = "cladogram")
axisPhylo(side=1)
dev.off()
system("open -a Skim.app ape.pdf")

pdf("ape.pdf", 15, 15)
plot.phylo(tree,use.edge.length = F,no.margin = T,type = "cladogram")
axisPhylo(side=1)
dev.off()
system("open -a Skim.app ape.pdf")




#require(popcorn)
source("scripts/popcorn_treemix.R")
tree <- read_treemix("tree.out.0")
pdf("pop.pdf", 15, 15)
p <- plot_treemix(tree) + theme_treemix()
print(p)
dev.off()
system("open -a Skim.app pop.pdf")


###########################################
######### Input(s) and Parameters #########
###########################################



###########################################
############### Wraggling #################
###########################################



###########################################
############### Analysis ##################
###########################################



###########################################
############### Plotting ##################
###########################################