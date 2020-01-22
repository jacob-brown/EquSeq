# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-21
# Last Modified: 2020-01-21

# Desc: read in each of the xml files and collate into one table, csv output

library(XML)
library(tidyverse)


xml_gather <- function(starting_level, id){

	# starting_level - of xml 
	# id - of BioProject

	# prepare the paths
	level = sprintf("//%s", starting_level) # level of xml to work on
	xml_doc <- sprintf("../data/esearch_xml/info_%s.xml", id) #filepath

	# parse and create dataframe
	doc <- xmlParse(xml_doc)
	xmldf <- xmlToDataFrame(nodes = getNodeSet(doc, level))

	return(xmldf)

}

a = xml_gather('DocumentSummary','PRJEB10098')

b = xml_gather('Id','PRJEB10098')


starting_level <- 'BioSample'
id <- 'PRJEB10098'

t <- read_xml(xml_doc)
top <- xml_find_all(t, "//BioSample")
top %>% 
    map(xml_attrs) %>% 
    map_df(~as.list(.)) %>%
    as.data.frame()







xmlToDataFrame(nodes = getNodeSet(doc, '//BioSample'))

myData = xmlToDataFrame(doc, stringsAsFactors = FALSE,) %>% 
  mutate_all(~type.convert(., as.is = T))


myXML = xmlParse(xml_doc)
myData = xmlToDataFrame(myXML, stringsAsFactors = FALSE,) %>% 
                        mutate_all(~type.convert(., as.is = T))

#####
library(xml2)
x <- read_xml(xml_doc)
x

xml_name(x)
xml_children(x)
xml_text(x)
a = xml_find_all(x, ".//SampleData")

a





#####


doc <- xmlParse(xml_doc)

bind_rows(xpathApply(doc, "//DocumentSummary", function(x) {
  parent <- data.frame(as.list(xmlAttrs(x)), stringsAsFactors=FALSE)
  kids <- bind_rows(lapply(xmlChildren(x), function(x) as.list(xmlAttrs(x))))
  cbind.data.frame(parent, kids, stringsAsFactors=FALSE)
}))


