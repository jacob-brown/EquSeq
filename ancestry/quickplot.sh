#! /bin/bash
while getopts "pakt" opt; do
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
    t) # plot 0 migratin tree and res
        Rscript ancestry/trees_plot.R s 
        open -a Skim.app results/ancestry/tree.pdf
        ;;
    \?) echo "usage [-p][-a]"
      ;;
  esac
done