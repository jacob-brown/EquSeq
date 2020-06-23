# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-09
# Last Modified: 2020-06-23

# Desc: Validate K for admixture analysis

###########################################
######### Input(s) and Parameters #########
###########################################

#location <- "results/ancestry/eu_wg_5kb_02maf/"
location <- "results/ancestry/eu_more_wg_5kb_05maf/"
command <- paste("sh ancestry/validate_K.sh ", location)
system(command)
logs <- as.data.frame(read.table(paste0(location, "/logfile.merge")))
logs <- abs(logs)
logs <- logs[order(logs$V1),]


###########################################
############### Plotting ##################
###########################################

pdf(file='results/ancestry/k_validation.pdf', 6, 6)
plot(logs$V1, logs$V2, type = "b", frame = FALSE, pch = 19, xlab = "K", ylab = "liklihood")
invisible(dev.off())
system("open -a Skim.app results/ancestry/k_validation.pdf")