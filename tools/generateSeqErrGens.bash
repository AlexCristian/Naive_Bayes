#! /bin/bash

# $1 is directory to operate on

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]
  then
    echo "usage: ./generateSeqErrGens.bash [source directory with base genomes] [destination directory] [with seq errs: true, false]
    WARNING: K-mer generator may not work for small spliced genome sizes.
    If this is the case, the '.kmr' files will be empty.
    "
    exit
fi

ORIGIN=$(dirname "$1")
TEST_GENS_FOLDER="$2"
N_FILES_EXPANDED_FROM_FAS=100

echo "Sequence all genomes in $1 [IN PROGRESS]"
echo "Data will be dumped to designated target folder $TEST_GENS_FOLDER"

if [ -d $TEST_GENS_FOLDER ] || [ -f $TEST_GENS_FOLDER ]
then
  rm -r $TEST_GENS_FOLDER
fi
mkdir $TEST_GENS_FOLDER

file_count=$(find $1 -iname "*.fasta" | wc -l)
progress=1
orig_dir=$(pwd)
files=$(find $1 -iname "*.fasta")
for file in $files
do
  echo "($progress/$file_count) Simulate read -  genome $file [IN PROGRESS]"
  filename="${file%.*}"
  class="$(dirname $filename)"
  class=${class##*/}
  filename="${filename##*/}"

  class_path=$TEST_GENS_FOLDER/$class
  if [ -f $class_path ]
  then
    rm $class_path
    mkdir $class_path
  fi
  if ! [ -d $class_path ]
  then
    mkdir $class_path
  fi

  if [ $3 == "true" ]
  then
    grinder -md poly4 3e-3 3.3e-8 -un 1 -rd 100 -tr 100 -rf $file -od $class_path > /dev/null
  else
    grinder -rd 100 -tr 100 -un 1 -rf $file -od $class_path > /dev/null
  fi
  mv $class_path/grinder-reads.fa $class_path/$filename.fasta
  rm $class_path/grinder-ranks.txt

  echo "($progress/$file_count) Simulate read - genome $file [COMPLETE]"
  progress=$(($progress+1))
done

progress=1

for file in $(find $TEST_GENS_FOLDER -iname "*.fasta")
do
  echo "($progress/$file_count) Expand FASTA file, $file [IN PROGRESS]"
  filename="${file%.*}"
  awk -v filename="$filename" -v reps=$N_FILES_EXPANDED_FROM_FAS 'BEGIN{OFS=""; ORS=""} />/{x=filename"."++i".fasta"; print $0"\n" > x; next;} {print tolower($1) > x; if(i == reps) exit;}' $filename.fasta
  rm $filename.fasta
  echo "($progress/$file_count) Expand FASTA file, $file [COMPLETE]"
  progress=$(($progress+1))
done

echo "Evolve all genomes in $1 [COMPLETE]"
