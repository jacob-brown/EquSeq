

cd dependancies
git clone https://github.com/samtools/htslib.git
git clone https://github.com/ANGSD/angsd.git 
cd htslib;make;cd ../angsd ;make HTSSRC=../htslib

cd ../
git clone https://github.com/Rosemeis/pcangsd.git pcangsd
python3 pcangsd/setup.py build_ext --inplace

# check if functional
# python pcangsd.py -h