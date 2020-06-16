#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-10
# Last Modified: 2020-06-10
# Desc: 


module load gsl

cd $EPHEMERAL/dependencies
git clone --recursive https://github.com/mfumagalli/ngsTools.git
cd ngsTools
export PKG_CONFIG_PATH=$GSL_HOME/lib/pkgconfig/
make


# local
cd ~/Desktop/
git clone --recursive https://github.com/mfumagalli/ngsTools.git
cd ngsTools
make

git clone https://github.com/fgvieira/ngsF-HMM.git
cd ngsF-HMM
make

	# if error? add the following
	export LIBRARY_PATH=/usr/local/Cellar/gsl/2.6/lib/
	export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:/usr/local/opt/openssl/lib/pkgconfig

	export LIBRARY_PATH="$LIBRARY_PATH:/usr/local/lib"