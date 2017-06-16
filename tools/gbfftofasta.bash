#! /bin/bash

function clearLastLine {
  tput cuu 1 && tput el
}

if [ -z "$1" ] || [ -z "$2" ]
  then
    echo "usage: ./gbfftofasta.bash [source directory with gbff files] [destination for FASTA files]
    "
    exit
fi

POOL_FOLDER=$2

if [ -d $POOL_FOLDER ] || [ -f $POOL_FOLDER ]
then
  rm -r $POOL_FOLDER
fi
mkdir $POOL_FOLDER

echo "Generate FASTA files [IN PROGRESS]"

progress=1
file_count=$(find $1 -iname "*.gbff" | wc -l)
for file in $(find $1 -iname "*.gbff")
do
  echo "($progress/$file_count) Generate FASTA for $file [IN PROGRESS]"
  filename="${file%.*}"
  filename="${filename##*/}"
  awk -v filename="$filename" ' BEGIN{TRUNC=0; OFS=""; ORS=""} /\/\//{if(TRUNC==1) {TRUNC=0; print "\n";}} TRUNC==1{$1=""; print $0} /^\s*ORIGIN\s*$/{TRUNC=1; print ">"filename"."++i"\n";}' $file > $POOL_FOLDER/$filename.fasta
  clearLastLine
  echo "($progress/$file_count) Generate FASTA for $file [COMPLETE]"
  progress=$(($progress+1))
done

clearLastLine
echo "Generate FASTA files [COMPLETE]"
