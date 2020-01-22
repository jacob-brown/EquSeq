#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-21
# Last Modified: 2020-01-21
# Desc: 

prj_modern=$( cat ../data/prj_modern.csv | head -5)
start="https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession="
end="&result=read_run"

for var in $prj_modern;
	do
		echo $var
		
		# save the esearch xml results ncbi
		esearch -db bioproject -query $var | elink -target biosample | efetch -format docsum > ../data/esearch_xml/info_$var.xml
		
		# save tab-delimenated results from ebi
		wget -O ../data/esearch_xml/$var.txt "${start}${var}${end}"

	done

exit

