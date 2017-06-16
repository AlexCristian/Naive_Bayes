#! /bin/bash

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]
  then
    echo "usage: ./copyFractionOfSeqs.bash [source directory with fasta files] [destination folder] [fraction 0-10]
    "
    exit
fi

echo "Copy files [IN PROGRESS]"

if ! [ -d $2 ]
then
  mkdir $2
fi

progress=1
file_count=$(find $1 -iname *.fasta | wc -l)
for file in $(find $1 -iname *.fasta)
do
  echo "($progress/$file_count) Iterating over files, $file [IN PROGRESS]"
  filename="${file%.*}"
  if [ $(( $RANDOM % 10 )) -gt $3 ]
  then
    class=$(dirname "$file")
    class="${class##*/}"
    filename="${file##*/}"
    if ! [ -d $2/$class ]
    then
      mkdir $2/$class
    fi
    cp $file $2/$class/$filename
  fi
  progress=$(($progress+1))
done
