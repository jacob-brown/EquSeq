#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-08
# Last Modified: 2020-04-08
# Desc: useful defind unix functions


###########################################
# timer function

time_start=$SECONDS

function timer {
	duration=$(($SECONDS - $time_start))
	echo -e "\n..........................\n"
 	echo "Time elapsed: " $duration " sec"
 	echo -e "\n..........................\n"
}

