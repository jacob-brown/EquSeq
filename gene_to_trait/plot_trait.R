# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-29
# Last Modified: 2020-07-15

# Desc: 


###########################################
################# Modules #################
###########################################

require(tidyverse)
#require(pdftools)

###########################################
############## Function(s) ################
###########################################



###########################################
######### Input(s) and Parameters #########
###########################################

df <- read.csv("results/gene_to_trait/genosub.csv")
counts_df <- read.csv("results/gene_to_trait/markercount.csv")
pop <- read.table("results/ancestry/clusters", header = TRUE, sep="\t")
pop$ind <- paste0("Ind", pop$IID-1)

cols <-  colnames(df)
col_inds <- cols[8:length(cols)]
col_info <- cols[2:7]
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
					"#0072B2", "#D55E00", "#CC79A7", "#000000")
borderCols <- c("white", "#56B4E9", "#009E73", "#E69F00", "#0072B2", 
					"#D55E00", "#CC79A7","#F0E442", "#000000")
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
								names_to = "ind", values_to = "GP")%>%
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
								names_to = "ind", values_to = "GP")%>%
					mutate(type = "homo_min", 
						ind = gsub(".2", "", ind, fixed = T))%>%
					select(-pres_homo_min)%>% 
					left_join(pop, by="ind")



#################
### homo major ###
homo_maj_df <- cbind(data, df[homo_maj]) %>% 
	select(-pres_homo_min, -pres_hetero)

# convert to NA if varient is absent

##### SHOULD GP FOR ABSENT ALLELES...
	# be assumed as 1 
	# or should they be noted as NA
	# or remain the same
#homo_maj_df[!homo_maj_df$pres_homo_maj,9:ncol(homo_maj_df)]# <- 0

homo_maj_df_piv <- tibble(homo_maj_df) %>%
					pivot_longer(cols = starts_with("Ind"), 
								names_to = "ind", values_to = "GP") %>%
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


###########################################
############### Plotting ##################
###########################################


#####################
### assign a type ###
	# hetero, homo min, homo maj


# create a group by of the max values at each marker-ind
filter_group <- data_all %>%
	group_by(clust_ind, phen_marker) %>%
	summarise(GP = max(GP)) %>%
	ungroup() %>%
	filter(GP != 0)


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


# check position of absent alleles and corret types
	# this will keep the GP, but highlight that the value it absent
if(!unique(homo_maj_df$absent == homo_min_df$absent & 
				hetero_df$absent == homo_min_df$absent)[1]){
	stop("absent values are not the same, issues with df positioning.")
}else{
	types[homo_maj_df$absent, ] <- "absent"

}


# value  matrix to long df
vdf <- as.data.frame(columnNameILike=row.names(vmat), vmat)
vdf$marker_phen <- rownames(vdf)

vdf_long <- tibble(vdf) %>%
				pivot_longer(cols = starts_with("Ind"), 
					names_to = "ind", values_to = "GP")

# category matrix to long df
tdf <- as.data.frame(columnNameILike=row.names(types), types)
tdf$marker_phen <- rownames(tdf)

tdf_long <- tibble(tdf) %>%
				pivot_longer(cols = starts_with("Ind"), 
					names_to = "ind", values_to = "type")


##### ABSENT SHOULD BE SWITCHED AROUND HERE ####
	# commented out below as 0 does no longer correspond to 0

# join  
tvdf <- left_join(tdf_long, vdf_long, by=c("marker_phen", "ind")) %>%
			left_join(pop, by="ind")# %>%
			#mutate(type = ifelse(GP == 0, "absent", as.character(type)))

tvdf$clust_ind <- paste0(tvdf$CLUSTER, "_", tvdf$ind)

# levels
clust_ind_levels <- sort(unique(tvdf$clust_ind))
clust_ind_levels <- rev(c("BENSON_Ind171", 
						clust_ind_levels[clust_ind_levels != "BENSON_Ind171"]))
tvdf$clust_ind <- factor(tvdf$clust_ind, levels=clust_ind_levels)

# make all min and maj labels homozygous
tvdf$type <- as.character(tvdf$type)
tvdf$type[tvdf$type %in% c("maj", "min", "min-maj")] <- "homo"
tvdf$type[tvdf$type %in% c("het-min", "het-maj")] <- "het-hom"  


#######################
##### Adjust DF #######
	# based on matteos feedback
		# equal GLs? no data --> NA
			# none of the data appears to be like this 
		# homo-het? 1 read --> NA
		# 0 or 1 read? --> NA
	# quick and dirty fix

# filter issue IDs
tvdf <- filter(tvdf, ind != "Ind29")

# 0-1 reads?
count_longdf <- tibble(counts_df) %>%
					dplyr::select(-X) %>%
					pivot_longer(cols = starts_with("Ind"), 
						names_to = "ind", values_to = "count") %>%
					mutate(marker = as.character(marker))

# equal GLs?  
	# none present
# homo-het?

# add NAs to:
	# unassigned genotypes
	# 0-1 reads from counts
	# 


tvdf$marker <- as.vector(
					sapply(tvdf$marker_phen, 
							function(x) str_split(x, "_", n=2)[[1]][2]))

alt_tvdf <- tvdf %>%
			mutate(type = ifelse(type %in% c("het-hom", "het-min-maj"), NA, type)) %>%
			left_join(count_longdf, by=c("marker", "ind")) %>%
			mutate(type = ifelse(count <= 1, NA, type),
					GP = ifelse(is.na(type), NA, GP)) %>%
			inner_join(df[c("phen","marker")], by="marker")
			#filter(!is.na(phen))
			# slight joining issue with na duplicates, corrected

###
# higher grouping (coat, disease, etc.)

	#write.table(unique(alt_tvdf$phen), "data/gene_variants/phen_raw.csv", row.names=F, quote=F,sep="\t") 
	# cp data/gene_variants/phen_raw.csv data/gene_variants/phen_adj.csv
		# manually assign higher categories

higher_grp <- read.table("data/gene_variants/phen_adj.csv",sep="\t", header=T)
colnames(higher_grp) <- c("phen", "phen_grp")

# final modifications before plot
	# add higher categories
	# change type to present or absent
plt_tvdf <- alt_tvdf %>%
				left_join(higher_grp, by="phen") %>% 
				mutate(marker_phen2 = paste0(phen, " (", marker, ")"),  
						type2  = 
							ifelse(type %in% c("het", "homo"), "present", type))

plt_tvdf[is.na(plt_tvdf$type2),]$type2 <- "ungenotyped"

#plt_tvdf <- filter(plt_tvdf, CLUSTER %in% as.character(head(sort(unique(tvdf$CLUSTER))))) 


plotGene <- function(){
# plot and save
g <- ggplot(plt_tvdf, mapping = aes(marker_phen2, ind, fill = GP)) + 
		facet_grid(vars(CLUSTER), vars(phen_grp), 
				scales = "free", space = "free", 
				switch="y")+
  		geom_tile()+
  		geom_point(aes(shape = type2), 
  				colour="white", size = 1.8, stroke=1)+
    	scale_fill_gradient(name  = "Probability", 
    				low="grey80", high="grey40", na.value = "#DC143C")+
    	scale_shape_manual(name  = "Classification", 
    			values=c(8, 4, NA), 
    			breaks=c("ungenotyped", "absent", "present"))+
    	ylab("")+
    	xlab("Phenotype")+
    	theme_classic()+
    	theme(axis.text.x = element_text(angle = 45, hjust=1, size=5),
    			strip.text.x = element_text(angle = 55, size = 7),
    			strip.text.y.left = element_text(angle = 0),
    			strip.placement.y = "outside",
    			strip.placement.x = "outside",
    			strip.background.x=element_blank(),
    			aspect.ratio = 1,
    			panel.spacing.x=unit(0.3, "lines") , 
    			panel.spacing.y=unit(0.1,"lines"),
    			legend.key=element_rect(fill="grey40")
    			)
    	#scale_fill_distiller(palette = "Blues", direction = 1, na.value =  "#DC143C")+
    	 
ggsave(filename="results/gene_to_trait/heat.comb.pdf", 
		plot=g, device ="pdf", width=10, height=30, units="in", dpi ="screen")
system("open -a Skim.app results/gene_to_trait/heat.comb.pdf")
}
plotGene()






#g <- ggplot(plt_tvdf, mapping = aes(marker_phen, clust_ind, fill = GP)) + 
#  		geom_tile(aes(colour=type, width=0.8, height=0.85), size=0.7) +
#    	scale_fill_gradient(low="white", high="blue", na.value = "red")+
#    	coord_fixed()+
#    	scale_colour_manual(values=c(borderCols), na.value = "red")+
#    	theme(axis.text.x = element_text(angle = 90, hjust=1, size=6))+
#    	guides(color=guide_legend(override.aes=list(fill=NA)))