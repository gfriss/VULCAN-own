# this script prepares the folder for the output
# runs the parallel python code
# them deletes the temporary folders used

# checking whether there has been a previous run and deletes that folder
if [ -d /tmp/datastore/s2555875/VULCAN_many ]
then
    rm -rf /tmp/datastore/s2555875/VULCAN_many
fi

mkdir /tmp/datastore/s2555875/VULCAN_many

cd /tmp/datastore/s2555875/VULCAN_many


for i in $(seq 0 $1); # only does until 99 now as a test
do
    if [[ $i < 10 ]]
    then    
        mkdir "sim_0${i}"
    else
        mkdir "sim_${i}"
    fi
done

cd /home/s2555875/
#python run_parallel.py