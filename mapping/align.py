#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-14
# Last Modified: 2020-04-14



"""  """

###########################################
################# Modules #################
###########################################

import argparse
import subprocess
import time

###########################################
################# Options #################
###########################################

parser = argparse.ArgumentParser(description='Align fastq to reference sequence.')
# input files
parser.add_argument("-i", "--input", dest="infile", nargs='+', type=str,
					required=True, help="read IN_FILE", metavar="IN_FILE")
# output name
parser.add_argument("-o", "--output", dest="outfile", type=str,
                  required=True, help="write OUT_FILE", metavar="OUT_FILE")
parser.add_argument("-d", "--directory", dest="dir", type=str,
                  required=True, help="file directory")
# reference genome
parser.add_argument("-r", "--ref", dest="refgen", type=str,
                  required=True, help="reference genome")
# threads
parser.add_argument("-t", "--threads", dest="threads", type=int,
                  required=False, help="number of threads")
# define args
args = parser.parse_args()

###########################################
############## Function(s) ################
###########################################

# timer
start = time.time() # start the timer from import

def timer():
	
	end = time.time()
	duration = end-start
	duration = round(duration)
	string = "\n..........................\n"\
			"   Time elapsed: {} sec"\
			"\n..........................\n"\
			.format(duration)

	print(string)


def bam(fastq1, fastq2, reference, outname, directory, threads=1):
	print('-----------------------')
	print("\nAlign sequences\n")
	# bwa mem $REF_GEN $FILE_1 $FILE_2 -t 31 > $DIR/aligned/$BASE_NAME'.sam'
	template = "bwa mem {ref} {file1} {file2} -t {nthreads}"\
				" > {dir}/aligned/{basename}.sam" 
	command = template.format(ref= reference, file1= fastq1,\
								 file2 = fastq2, nthreads= threads,\
								 dir=directory, basename=outname)
	#subprocess.run(command)
	print(command)


###########################################
######### Input(s) and Parameters #########
###########################################

fastq1 = args.infile[0]
fastq2 = args.infile[1]
args.outfile
args.refgen
args.dir
args.threads
subprocess.run(["echo", "hello", "&& echo"])


# bwa mem $REF_GEN $FILE_1 $FILE_2 -t 31 > $DIR/aligned/$BASE_NAME'.sam'
bam(args.infile[0], args.infile[1], args.refgen, args.outfile, args.dir, args.threads)

#bwa mem $REF_GEN $FILE_1 $FILE_2 -t 31 > $DIR/aligned/$BASE_NAME'.sam'



#echo '-----------------------'
#echo -e "\nConvert to bam\n"
#
#samtools view -bS --threads 31 $DIR/aligned/$BASE_NAME'.sam' > \
#		$DIR/converted/$BASE_NAME'.bam'
#
#
#
#echo '-----------------------'
#echo -e "\nSorting\n"
#
#samtools sort -m 60GiB --threads 31 $DIR/converted/$BASE_NAME'.bam' -o  \
#		$DIR/sorted/$BASE_NAME'.sorted.bam'
#
#
#
#echo '-----------------------'
#echo -e "\nIndex\n"
#
#samtools index $DIR/sorted/$BASE_NAME'.sorted.bam'
#
#
#
#echo '-----------------------'
#echo -e "\nFlagstat\n"
#
## should be high due to 
#samtools flagstat $DIR/sorted/$BASE_NAME'.sorted.bam' > \
#		$DIR/stats/$BASE_NAME'.stat.txt'



