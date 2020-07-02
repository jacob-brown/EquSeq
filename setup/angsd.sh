
# angsd
cd dependancies
git clone https://github.com/samtools/htslib.git
git clone https://github.com/ANGSD/angsd.git 
cd htslib;make;cd ../angsd ;make HTSSRC=../htslib

# check if functional
#angsd/angsd

# pcangsd
cd ../
git clone https://github.com/Rosemeis/pcangsd.git pcangsd
python pcangsd/setup.py build_ext --inplace

# check if functional
# python pcangsd.py -h

module load anaconda3/personal
python3 setup.py build_ext --inplace

# issues? try in conda
source activate myenv
python pcangsd.py -h
conda deactivate