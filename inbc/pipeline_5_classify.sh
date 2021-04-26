#!/bin/bash
### tell SGE to use bash for this script
#$ -S /bin/bash
### execute the job from the current working directory, i.e. the directory in which the qsub command is given
#$ -cwd
### join both stdout and stderr into the same file
#$ -j y
### set email address for sending job status
#$ -M zz374@drexel.edu
### project - basically, your research group name with "Grp" replaced by "Prj"
#$ -P rosenPrj
### select parallel environment, and number of job slots
#$ -pe shm 5
### request 15 min of wall clock time "h_rt" = "hard real time" (format is HH:MM:SS, or integer seconds)
#$ -l h_rt=48:00:00
### a hard limit 16 GB of memory per slot - if the job grows beyond this, the job is killed
#$ -l h_vmem=3.5G
### want nodes with at least 15 GB of free memory per slot
#$ -l m_mem_free=3G
### select the queue all.q, using hostgroup @intelhosts
#$ -q all.q
##$ -t 1-100:1

. /etc/profile.d/modules.sh
. ~/.bashrc
### These four modules must ALWAYS be loaded
module load shared
module load proteus
module load sge/univa
module load gcc
module load boost/openmpi/gcc/64/1.57.0

which python

KSIZE=$(cat kmer.txt)

FILE=airwaystest/short_read/2017.12.04_18.56.22_sample_4/reads/anonymous_reads.fq.gz
total=$(zcat $FILE | wc -l)
index=$SGE_TASK_ID

delta=10000
start=$(((index - 1)*delta*4 + 1)) #40001
end=$((start + delta*4 - 1))
next=$((end + 1))
if [ $end -gt $total ]
then
end=$total
next=$total
fi
site_dir=testing-$(echo $FILE | cut -d/ -f1)
sample_dir=$(echo $FILE | cut -d/ -f3)-fold-$start
output_dir=$site_dir/$sample_dir
## mkdir $site_dir
mkdir $output_dir
output_file=$output_dir/tmp.txt

echo "${start},${end}p;${next}q"
echo $output_file

zcat $FILE | sed -n "${start},${end}p;$((end + 1))q" > $output_file
bash 2_fq2fa.sh $output_file

mkdir $output_dir/read_file
which python
echo ${output_file%.txt}.fna
echo $output_dir/read_file
SRC=$output_dir

python split_reads.py ${output_file%.txt}.fna $output_dir/read_file

rm $output_dir/tmp.*

module load jellyfish/gcc/2.2.10

start=`date +%s`
bash jellyfish_gen_new.bash $SRC $KSIZE false
end=`date +%s`
runtime=$((end-start))
echo "jellyfish runtime is $runtime"

CORE_NUM=5
CODE=$(pwd)/nb-classify-batch.bash
MODEL=/lustre/scratch/zz374/model_$KSIZE
start=`date +%s`
$CODE $SRC $MODEL $KSIZE $((CORE_NUM - 1))
end=`date +%s`
runtime=$((end-start))
echo "classification time is $runtime"
rm -rf $SRC
