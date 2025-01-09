#!/bin/bash

# pe request
#$ -pe mpi_1 1

# our Job name

#$ -S /bin/bash

#$ -q all.q

#$ -V
#$ -cwd


export OMP_NUM_THREADS=$NSLOTS
WDIR=/work/hkim/gaussian/compound_2_ketone_1/t1_opt/B3LYP/UB3LYP/duschinsky_IC/rate.$$
DIR=`pwd`
MOL='duschinsky_read'
mkdir -p $WDIR
cp $MOL.py $WDIR
cp *.txt $WDIR
cd $WDIR

start_time=$(date +"%Y-%m-%d %H:%M:%S")
python $MOL.py > $MOL_second_half.out
cp * $DIR
rm -rf $WDIR
end_time=$(date +"%Y-%m-%d %H:%M:%S")
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
elapsed_seconds=$((end_seconds - start_seconds))

echo "Job started at: $start_time"
echo "Job ended at: $end_time"
echo "Elapsed time: $elapsed_seconds seconds"
echo "Used $NSLOTS cores"

