# Usage: Rscript -i infile.covar -c component1-component2 -a annotation.file -o outfile.pdf

library(methods)
library(optparse)
library(ggplot2)

# colorblind friendly: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# To use for fills, add: scale_fill_manual(values=cbPalette)
# To use for line and point colors, add: scale_colour_manual(values=cbPalette)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file (output from ngsCovar)'),
                    make_option(c('-c','--comp'), action='store', type='character', default=1-2, help='Components to plot'),
                    make_option(c('-a','--annot_file'), action='store', type='character', default=NULL, help='Annotation file with individual classification (2 column TSV with ID and ANNOTATION)'),
                    make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file')
                    )
opt <- parse_args(OptionParser(option_list = option_list))

# Annotation file is in plink cluster format

#################################################################################

# Read input file
covar <- read.table(opt$in_file, stringsAsFact=F);

# Read annot file
annot <- read.table(opt$annot_file, sep="\t", header=T); # note that plink cluster files are usually tab-separated

# Parse components to analyze
comp <- as.numeric(strsplit(opt$comp, "-", fixed=TRUE)[[1]])

# Eigenvalues
eig <- eigen(covar, symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val, digits=3)*100,"\n");

# Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))
PC$Pop <- factor(annot$CLUSTER)

title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="",collapse="")

x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")

ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop")) + ggtitle(title) + scale_colour_manual(values=cbPalette)
ggsave(opt$out_file)
unlink("Rplots.pdf", force=TRUE)