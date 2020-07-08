#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-16
# Last Modified: 2020-07-07
# Desc: remove snps from treemix file that have 0,0 for novel sample
	# treemix returns errors whentoo many are present in one sample


import subprocess, gzip

p = subprocess.Popen(["zcat < treemix.frq.gz"], stdout=subprocess.PIPE, \
		stderr=subprocess.PIPE, shell = True)
out, er = p.communicate()
out_decode = out.decode()
pos_str = out_decode.split("\n")
pos_str.remove("")
data = [i.split(" ") for i in pos_str]
bensonindex = data[0].index("BENSON") # sample name
dataout = [" ".join(i) for i in data if i[bensonindex] != "0,0"] # sites to remove
content = "\n".join(dataout).encode()
with gzip.open("treemix.benson.frq.gz", "wb") as f:
    f.write(content)