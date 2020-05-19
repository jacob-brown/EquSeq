#! /bin/bash
while getopts "pa" opt; do
  case ${opt} in
  	p) # pca
		Rscript ancestry/plotPCA.R
		open -a Skim.app results/ancestry/ALL.PCA.pdf
		;;
    a) # admixture
       	Rscript ancestry/plotAdmix.R 
       	open -a Skim.app results/ancestry/ALL.MIX.pdf
      	;;
    \?) echo "usage [-p][-a]"
      ;;
  esac
done