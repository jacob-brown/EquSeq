# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-26
# Last Modified: 2020-06-26

# Desc: 


###########################################
################# Modules #################
###########################################


###########################################
############## Function(s) ################
###########################################



###########################################
######### Input(s) and Parameters #########
###########################################

files <- list.files()
files_to_use <- files[grep(".sfs", files)]

###########################################
############### Plotting ##################
###########################################

len <- length(files_to_use)

pdf("sfs.pdf", 10, 5)
par(mfrow=c(2,len/2))
for(elem in files_to_use){
	sfs<-scan(elem)
	barplot(sfs[-1], main=elem)
}
dev.off()
system("open -a Skim.app sfs.pdf")


