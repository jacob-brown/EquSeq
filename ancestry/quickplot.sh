#! /bin/bash
while getopts "pakt" opt; do
  case ${opt} in
  	p) # pca
		Rscript ancestry/plotPCA.R
		;;
    a) # admixture
       	Rscript ancestry/plotAdmixPopHelper.R 
       	#open -a Skim.app results/ancestry/ALL.MIX.pdf
      	;;
    k) # validate K
        Rscript ancestry/validate_K.R 
        open -a Skim.app results/ancestry/k_validation.pdf
        ;;
    t) # plot 0 migratin tree and res
        Rscript ancestry/plot_trees.R s 
        open -a Skim.app results/ancestry/tree.pdf
        ;;
    \?) echo "usage [-p][-a]"
      ;;
  esac
done