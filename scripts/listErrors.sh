# usage
	# list errors in file with specific grep string
	# sh listErrors.sh fastqDownload.sh closed

FILENAME=$1

files=($(ls $FILENAME.e*)) 
errorGrep="$2"

echo "printing problem files with grep: " $errorGrep

for i in "${files[@]}"
do
	if grep -q $errorGrep $i; then
	    echo $i
	fi
done


