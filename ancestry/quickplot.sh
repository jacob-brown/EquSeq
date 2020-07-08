#! /bin/bash
while getopts "pak" opt; do
  case ${opt} in
  	p) # pca
		Rscript ancestry/plotPCA.R
		;;
    a) # admixture
       	Rscript ancestry/plotAdmix.R 
       	open -a Skim.app results/ancestry/ALL.MIX.pdf
      	;;
    k) # validate K
        Rscript ancestry/validate_K.R 
        open -a Skim.app results/ancestry/k_validation.pdf
        ;;
    \?) echo "usage [-p][-a]"
      ;;
  esac
done