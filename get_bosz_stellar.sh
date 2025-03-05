# script downloading the necessary synthetic spectra from the BOSZ library

# declare arrays and values for Teff, logg
declare -a Teff=("t2800" "t3100" "t3400" "t3700" "t4000" "t4250" "t4500" "t4750" "t5000" "t5250" "t5500" "t5750" "t6000" "t6250" "t6500")
declare -a logg=("g+5.0" "g+5.0" "g+5.0" "g+5.0" "g+4.5" "g+4.5" "g+4.5" "g+4.5" "g+4.5" "g+4.5" "g+4.5" "g+4.5" "g+4.0" "g+4.0" "g+4.0")

# get length of an array
arraylength=${#Teff[@]}
# change directory to the scratch directory
cd /scratch/s2555875/BOSZ_spectra
# use for loop to download and extract the tar files of spectra
for (( i=0; i<${arraylength}; i++ ));
do
    wget https://archive.stsci.edu/hlsps/bosz/bosz2024/r500/m+0.00/bosz2024_mp_${Teff[$i]}_${logg[$i]}_m+0.00_a+0.00_c+0.00_v0_r500_resam.txt.gz
    gunzip bosz2024_mp_${Teff[$i]}_${logg[$i]}_m+0.00_a+0.00_c+0.00_v0_r500_resam.txt.gz
done
cd /home/s2555875/VULCAN-2