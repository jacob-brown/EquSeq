#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-12
# Last Modified: 2020-08-20
# Desc: check log files for errors


# zip files
	# zip -r fastqDown.zip fastqDownload.sh.*

mkdir sandbox/fasta_logs

cd sandbox/fasta_logs
mkdir log_pass
mkdir log_fail

scp ADD_YOUR_SCP_HERE .

unzip fastqDown.zip

files=($(ls fastqDownload.sh.e*)) 
errorGrep="curl"

echo "assigning to pass/fail"

for i in "${files[@]}"
do
	outfile=($(echo $i | sed 's/.e/.o/g'))
	if grep -q $errorGrep $i; then
	    mv $i $outfile log_fail/
	else
		mv $i $outfile log_pass/
	fi
done

echo "number of files failed"
ls log_fail/*.e* | wc -l

echo "getting the error RUN codes"
cd log_fail

fail=($(ls fastqDownload.sh.o*)) 


(for f in "${fail[@]}";do awk '/^requesting fasta: / {print $3}' $f; done) > issue.run.codes

cat issue.run.codes

### generate an info table for the error codes 
python3 scripts/fastq_download_check.py

### rerun fastqDownload script - see comments ###
