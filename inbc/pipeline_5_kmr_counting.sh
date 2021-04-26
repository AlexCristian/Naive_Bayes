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
#$ -pe shm 30
### request 15 min of wall clock time "h_rt" = "hard real time" (format is HH:MM:SS, or integer seconds)
#$ -l h_rt=120:00:00
### a hard limit 8 GB of memory per slot - if the job grows beyond this, the job is killed
#$ -l h_vmem=2G
### want nodes with at least 6 GB of free memory per slot
#$ -l m_mem_free=1.5G
### select the queue all.q, using hostgroup @intelhosts
#$ -q long.q
#$ -t 1-7:1

. /etc/profile.d/modules.sh

### These four modules must ALWAYS be loaded
module load shared
module load proteus
module load sge/univa
module load gcc

module load jellyfish/gcc/2.2.10

kmer=$(cat kmer.txt)
DIGIT=9
TENTH=$(($SGE_TASK_ID - 1))
SRC=test_set/test_fold_$TENTH$DIGIT/test
echo $SRC
start=`date +%s`
bash jellyfish_gen.bash $SRC $kmer false
end=`date +%s`
runtime=$((end-start))
echo "jellyfish runtime is $runtime"
