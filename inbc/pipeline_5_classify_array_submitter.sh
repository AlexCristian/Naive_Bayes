#!/bin/bash
. /etc/profile.d/modules.sh

module load shared
module load proteus
module load sge/univa

# 1:3333

FILE=airwaystest/short_read/2017.12.04_18.56.22_sample_4/reads/anonymous_reads.fq.gz
total=$(zcat $FILE | wc -l)

site_dir=testing-$(echo $FILE | cut -d/ -f1)

mkdir $site_dir

##w1=$( qsub -terse -t 1:100:1 pipeline_5_classify.sh | cut -f1 -d. )
##w2=$( qsub -terse -t 101:200:1 -hold_jid $w1 pipeline_5_classify.sh | cut -f1 -d.)
##w3=$( qsub -terse -t 201:300:1 -hold_jid $w2 pipeline_5_classify.sh | cut -f1 -d.)
##w4=$( qsub -terse -t 301:400:1 -hold_jid $w3 pipeline_5_classify.sh | cut -f1 -d.)
##qsub -terse -t 401:500:1 -hold_jid $w4 pipeline_5_classify.sh

wid=$( qsub -terse -t 500:1000:1 pipeline_5_classify.sh | cut -f1 -d. )
echo $wid
wid=$( qsub -terse -t 1001:1500:1 -hold_jid $wid pipeline_5_classify.sh | cut -f1 -d.)
echo $wid
wid=$( qsub -terse -t 1501:2000:1 -hold_jid $wid pipeline_5_classify.sh | cut -f1 -d.)
echo $wid
wid=$( qsub -terse -t 2001:2500:1 -hold_jid $wid pipeline_5_classify.sh | cut -f1 -d.)
echo $wid
wid=$( qsub -terse -t 2501:3000:1 -hold_jid $wid pipeline_5_classify.sh | cut -f1 -d.)
echo $wid
wid=$( qsub -terse -t 3001:3333:1 -hold_jid $wid pipeline_5_classify.sh | cut -f1 -d.)
echo $wid
