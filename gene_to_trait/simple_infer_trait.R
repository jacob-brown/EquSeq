# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-29
# Last Modified: 2020-06-30

# Desc: 


###########################################
################# Modules #################
###########################################

require(tidyverse)
require(pdftools)

###########################################
############## Function(s) ################
###########################################



###########################################
######### Input(s) and Parameters #########
###########################################

df = read.csv("results/gene_to_trait/beaglesub.csv")
pop <- read.table("results/ancestry/clusters", header = TRUE, sep="\t")
pop$ind <- paste0("Ind", seq(1:nrow(pop))-1)

cols <-  colnames(df)
col_inds <- cols[8:length(cols)]
col_info <- cols[2:7]
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

###########################################
############### Wraggling #################
###########################################

# get columns
hetero <- col_inds[str_sub(col_inds, start= -2) == ".1"]
homo_min <- col_inds[str_sub(col_inds, start= -2) == ".2"]
homo_maj <- col_inds[!(col_inds %in% homo_min | col_inds %in% hetero)]


# prep the general data columns
df$new <- as.character(df$new)
df$allele1 <- as.character(df$allele1)
df$allele2 <- as.character(df$allele2)
df$marker <- as.character(df$marker)
df$phen <- as.character(df$phen)

data <- tibble(df[col_info]) %>%
			mutate(major = allele1,
					minor = allele2,
					absent = ifelse(new == major, FALSE,
										ifelse(new == minor, FALSE, TRUE))) %>%
			select(-ref, -allele1, -allele2) %>%
			mutate(pres_hetero = ifelse(!absent, TRUE, FALSE),
					pres_homo_min = ifelse(!absent & new == minor, TRUE, FALSE),
					pres_homo_maj = ifelse(!absent & new == major, TRUE, FALSE),
					phen_marker = paste0(phen, "_", marker))




# seperate data

##############
### hetero ###
hetero_df <- cbind(data, df[hetero]) %>% 
		select(-pres_homo_min, -pres_homo_maj)

# convert to NA if varient is absent
hetero_df[!hetero_df$pres_hetero,9:ncol(hetero_df)] <- 0




# heatmap
hetero_df_piv <- tibble(hetero_df) %>%
					pivot_longer(cols = starts_with("Ind"), 
								names_to = "ind", values_to = "GL")%>%
					mutate(type = "hetero", 
							ind = gsub(".1", "", ind, fixed = T))%>%
					select(-pres_hetero) %>% 
					left_join(pop, by="ind")


#################
### homo minor ###
homo_min_df <- cbind(data, df[homo_min]) %>% 
	select(-pres_hetero, -pres_homo_maj) 

# convert to NA if varient is absent
homo_min_df[!homo_min_df$pres_homo_min,9:ncol(homo_min_df)] <- 0

# ggplot heatmap
homo_min_df_piv <- tibble(homo_min_df) %>%
					pivot_longer(cols = starts_with("Ind"), 
								names_to = "ind", values_to = "GL")%>%
					mutate(type = "homo_min", 
						ind = gsub(".2", "", ind, fixed = T))%>%
					select(-pres_homo_min)%>% 
					left_join(pop, by="ind")




#################
### homo major ###
homo_maj_df <- cbind(data, df[homo_maj]) %>% 
	select(-pres_homo_min, -pres_hetero)

# convert to NA if varient is absent
homo_maj_df[!homo_maj_df$pres_homo_maj,9:ncol(homo_maj_df)] <- 0



# ggplot heatmap
homo_maj_df_piv <- tibble(homo_maj_df) %>%
					pivot_longer(cols = starts_with("Ind"), 
								names_to = "ind", values_to = "GL") %>%
					mutate(type = "homo_maj") %>%
					select(-pres_homo_maj) %>% 
					left_join(pop, by="ind")





####### 
### Combine datasets ###
	# and add clusters 
data_all <- rbind(homo_maj_df_piv, homo_min_df_piv, hetero_df_piv)
data_all$clust_ind <- paste0(data_all$CLUSTER, "_", data_all$ind)
###########################################
############### Analysis ##################
###########################################

#inds <- c(paste0("Ind", seq(0:44)), "ind177")
#tmp_data_all <- data_all %>% filter(ind %in% inds)

###########################################
############### Plotting ##################
###########################################


p <- ggplot(data_all, aes(clust_ind, phen_marker, fill= GL)) + 
 	 facet_grid(rows=vars(type))+
 	 geom_tile()+
 	 scale_fill_gradient(low="white", high="blue") +
 	 theme(axis.text.x = element_text(angle = 90))

#pdf("sandbox/heat.pdf", 20, 20)
#print(p)
#dev.off()
#system("open -a Skim.app sandbox/heat.pdf")

#ggsave(filename="results/gene_to_trait/heat.pdf", 
#		plot=p, device ="pdf", width=20, height=20, units="in", dpi ="screen")
#
#bitmap <- pdf_render_page("results/gene_to_trait/heat.pdf", page = 1)

# issues with R and ariel font with newest mac update... 
#png::writePNG(bitmap, "results/gene_to_trait/heat.png")
#system("open -a Preview.app results/gene_to_trait/heat.png")


#pdf("sandbox/heat2.pdf", 20, 7)
#
#for(i in 1:3){
#	df <- list(data_all, homo_min_df_piv, hetero_df_piv)[[i]]
#	type <- unique(df$type)
#
#	df$clust_ind <- paste0(df$CLUSTER, "_", df$ind)
#
#	g <- ggplot(df, aes(clust_ind, phen_marker, fill= GL)) + 
# 	 geom_tile()+
# 	 ggtitle(type)+
# 	 scale_fill_gradient(low="white", high="blue") +
# 	 theme(axis.text.x = element_text(angle = 90))
#
# 	print(g)
#}
#
#dev.off()
#system("open -a Skim.app sandbox/heat2.pdf")


#####################
### assign a type ###
	# hetero, homo min, homo maj


# create a group by of the max values at each marker-ind
filter_group <- data_all %>%
	group_by(clust_ind, phen_marker) %>%
	summarise(GL = max(GL)) %>%
	ungroup() %>%
	filter(GL != 0)
	#data.frame() %>%
	#tibble()
#	inner_join(data_all, by=c("clust_ind", "phen_marker", "GL")) 
#data_all %>%
#	mutate(type = ifelse(GL==0.333333, "ungeno", type)) %>% 
#	filter(type=="ungeno") %>% arrange(ind, phen_marker)
#	unique() %>%
#	inner_join(filter_group, 
#		by=c("clust_ind", "phen_marker", "GL"))



### convert to matrix ###

# hetero
hetero_mat <- as.matrix(hetero_df[9:ncol(hetero_df)])
rownames(hetero_mat) <- hetero_df$phen_marker
colnames(hetero_mat) <- gsub(".1", "", colnames(hetero_mat), fixed = T)

# homo minor
homo_min_mat <- as.matrix(homo_min_df[9:ncol(homo_min_df)])
rownames(homo_min_mat) <- homo_min_df$phen_marker
colnames(homo_min_mat) <- gsub(".2", "", colnames(homo_min_mat), fixed = T)

# homo major
homo_maj_mat <- as.matrix(homo_maj_df[9:ncol(homo_maj_df)])
rownames(homo_maj_mat) <- homo_maj_df$phen_marker

# compare
winners <- pmax(hetero_mat, homo_min_mat, homo_maj_mat)

# identify the winners
het_logic <- hetero_mat == winners
min_logic <- homo_min_mat == winners
maj_logic <- homo_maj_mat == winners

# determine which values are draws
draws <- ifelse(maj_logic + het_logic == 2, "het-maj",
			ifelse(maj_logic + min_logic == 2, "min-maj",
				ifelse(het_logic + min_logic == 2, "het-min",
					ifelse(het_logic + maj_logic +  min_logic == 2, 
							"het-min-maj", NA))))
# determine clear winners
types <- ifelse(!is.na(draws), draws,
			ifelse(het_logic, "het",
				ifelse(min_logic, "min",
					ifelse(maj_logic, "maj", NA))))


# condition matrix
maj_cond <- (types == "maj" | 
				types == "het-maj" | 
				types == "min-maj")

min_cond <- (types == "min" | 
				types == "het-min")

het_cond <- (types == "het")

# value matrix
vmat <- ifelse(maj_cond, homo_maj_mat, 
			ifelse(min_cond, homo_min_mat,
	 			ifelse(het_cond, hetero_mat,
					NA)))


# value  matrix to long df
vdf <- as.data.frame(columnNameILike=row.names(vmat), vmat)
vdf$marker_phen <- rownames(vdf)

vdf_long <- tibble(vdf) %>%
				pivot_longer(cols = starts_with("Ind"), 
					names_to = "ind", values_to = "GL")

# category matrix to long df
tdf <- as.data.frame(columnNameILike=row.names(types), types)
tdf$marker_phen <- rownames(tdf)

tdf_long <- tibble(tdf) %>%
				pivot_longer(cols = starts_with("Ind"), 
					names_to = "ind", values_to = "type")


# join  
tvdf <- left_join(tdf_long, vdf_long, by=c("marker_phen", "ind")) %>%
			left_join(pop, by="ind") %>%
			mutate(type = ifelse(GL == 0, "absent", as.character(type)))

tvdf$clust_ind <- paste0(tvdf$CLUSTER, "_", tvdf$ind)

#tvdf <- filter(tvdf, CLUSTER %in% c("BENSON", "Friesian dwarf", "Standardbred", "Akhal-Teke"))
#filter(tvdf, CLUSTER == "BENSON") %>% data.frame()
# plot and save
g <- ggplot(tvdf, mapping = aes(marker_phen, clust_ind, fill = GL)) + 
  		geom_tile(aes(colour=type, width=0.8, height=0.85), size=0.7) +
    	scale_fill_gradient(low="white", high="blue")+
    	coord_fixed()+
    	scale_colour_manual(values=c("white",cbbPalette))+
    	theme(axis.text.x = element_text(angle = 90, hjust=1, size=6))+
    	guides(color=guide_legend(override.aes=list(fill=NA)))
    	
    	 
ggsave(filename="results/gene_to_trait/heat.comb.pdf", 
		plot=g, device ="pdf", width=10, height=30, units="in", dpi ="screen")
system("open -a Skim.app results/gene_to_trait/heat.comb.pdf")
# issues with R and ariel font with newest mac update... 
#bitmap <- pdf_render_page("results/gene_to_trait/heat.comb.pdf", page = 1)
#png::writePNG(bitmap, "results/gene_to_trait/heat.comb.png")
#system("open -a Preview.app results/gene_to_trait/heat.comb.png")


