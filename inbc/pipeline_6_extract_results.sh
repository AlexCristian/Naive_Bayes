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
#$ -pe shm 1
### request 15 min of wall clock time "h_rt" = "hard real time" (format is HH:MM:SS, or integer seconds)
#$ -l h_rt=48:00:00
### a hard limit 8 GB of memory per slot - if the job grows beyond this, the job is killed
#$ -l h_vmem=3G
### want nodes with at least 6 GB of free memory per slot
#$ -l m_mem_free=2G
### select the queue all.q, using hostgroup @intelhosts
#$ -q all.q

. /etc/profile.d/modules.sh
. ~/.bashrc
### These four modules must ALWAYS be loaded
module load shared
module load proteus
module load sge/univa
module load gcc

echo -e '@@SEQUENCEID\tTAXID' > Outputfile.txt
for file in $(find results -iname "pipeline_5_classify.sh.o*")
do
file_name=tmp-$(echo $file | cut -d. -f4).txt
grep "Genome" $file | cut -d, -f3,2 > $file_name 
echo $file_name
python output_formatting.py $file_name
done
