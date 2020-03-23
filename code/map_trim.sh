#!/bin/bash
# trim fastq files

# catch input files 
read FILE_IN FILE_OUT
echo -e 'trimming: ' $FILE_IN '\nwith output:\n' $FILE_OUT 

# trim
echo '=================================='
echo -e "\nTrimming\n"


fastx_trimmer -l 90 -i $FILE_IN -o $FILE_OUT
		# ~ 8 mins per read


exit