
# Desc: Get files from ncbi sra


#----- load modules ----#
echo '=================================='
echo -e "\nLoad modules\n"
module load sra-toolkit/2.8.1

# run code
prefetch ERR466186
prefetch ERR868004

fastq-dump -X 5 -Z ERR466186





