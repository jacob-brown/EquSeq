FILENAME=$1

files=($(ls $FILENAME.sh.e*)) 


for i in "${files[@]}"
do
	tail -1 $i
	echo $i
done


# bash listErrors.sh sra_mapping.sh