#! /bin/bash
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=1:mem=1gb
# Desc: generate pass/fail statistic files from kraken output

# qsub -J 0-23 kraken_stat.sh
	# qsub -J 0-3 kraken_stat.sh

# 3 sec for 10000 lines
# each file has ~ 30000000 lines
# 1 file should be finished in ~2.5 hours


source $HOME/genomics/EquSeq/scripts/unix_functions.sh
module load anaconda3/personal
DIR=$EPHEMERAL/oral_diversity/
KFILES=$DIR/kraken_out/
ALL_FILES=($(ls $KFILES/*.kraken))

FILE=${ALL_FILES[$PBS_ARRAY_INDEX]} 
#FILE=$DIR/all.tmp.kraken

echo '=================================='
echo -e "\nreading file: " $FILE "\n"


python $HOME/genomics/EquSeq/oral_diversity/kraken_explore.py\
		-i $FILE \
		-o $DIR/stat.$PBS_ARRAY_INDEX \
		-t 0.2

echo -e "\nfile out: " $DIR/stat.$PBS_ARRAY_INDEX












### test speed ###
#head -n10000 kraken_out/V300044309_L2_B5GHORlfyRAAAAAAA-517.kraken > all.tmp.kraken
#wc -l all.tmp.kraken
#FILE=all.tmp.kraken

#time python $HOME/genomics/EquSeq/oral_diversity/kraken_explore.py -i $FILE -o $DIR/stat.$FILE -t 0.2

# 3 sec for 10000 lines
# each file has ~ 30000000 lines
# 1 file should be finished in ~2.5 hours

#COUNTER=0
#for file in $KFILES/*.kraken
#do
#			
#	python $HOME/genomics/EquSeq/oral_diversity/kraken_explore.py -i $file -o $DIR/stat.$COUNTER -t 0.2
#	
#	let COUNTER+=1
#	
#done


##file=all.tmp.kraken
#file=all.517.kraken
#WD=$(wc -l $file | awk {'print $1'}) # word count
#STEP=$(( $WD / 50 )) # steps
#STEPS=($(seq 0 $STEP $WD))
#
#echo '=================================='
#echo -e "\ncreating tmp files\n"
#
## create tmp files
#for i in {1..4}
#do
#	start="${STEPS[$i-1]}"
#	start_cor=$(( $start + 1 ))
#	end="${STEPS[$i]}"
#
#	echo -e "\nstart: " $start_cor " end: " $end 
#
#	echo -e "writing " tmp_files/tmp.$i.$file
#	sed ''"$start_cor"','"$end"'!d' $file > tmp_files/tmp.$i.$file
#
#done
#
#echo '=================================='
#echo -e "\nGenrating stat files\n"
#
## generate stat files
#for i in {1..4}
#do
#	(tmpfile=tmp_files/tmp.$i.all.517.kraken;
#	echo $tmpfile;
#	python $HOME/genomics/EquSeq/oral_diversity/kraken_explore.py -i $tmpfile  -o stat.$i -t 0.2) &
#done; wait
#timer


#wc < tmp_files/tmp.1.all.tmp.kraken | awk {'print $1'}
#echo '=================================='
#echo -e "\nclear tmp files\n"
#
#rm tmp_files/tmp.*.$file

