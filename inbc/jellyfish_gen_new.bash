#! /bin/bash

function clearLastLine {
  tput cuu 1 && tput el
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]
  then
    echo "usage: ./jellyfish_gen.bash [source directory with fna files] [k-mer size] [is full genome: true, false] [resume from #]
    "
    exit
fi

# if [ -z "$4" ]
# then
  seq_err=0
# else
#   seq_err=$4
# fi

echo "Generate kmer count [IN PROGRESS]"

# Call jellyfish
progress=1
file_count=$(find $1 -iname "*.fna" | wc -l)
for file in $(find $1 -iname "*.fna")
do
  if ! [ -z "$4" ] && [ $progress -lt $4 ]
  then
    progress=$(($progress+1))
    continue
  fi
  echo "($progress/$file_count) Generate kmer count for $file [IN PROGRESS]"
  filename="${file%.*}"
  size=$( cat "$file" | wc -c )
  size=$( awk -v a="$size" -v b="$seq_err" 'BEGIN{printf("%d\n",(a + a * b) + 0.5)}' )
  if [ $3 == "true" ]
  then
    sgn=1
    while [ $sgn != "0" ]
    do
      timeout 300 jellyfish count -m $2 -s $size -t 32 -o $filename.jf $file
      sgn=$?
      if [ $sgn != "0" ]
      then
        echo "Stuck. Retrying"
        rm $filename.jf
      fi
    done
  else
    sgn=1
    while [ $sgn != "0" ]
    do
      timeout 300 jellyfish count -m $2 -s $size -t 32 -C -o $filename.jf $file
      sgn=$?
      if [ $sgn != "0" ]
      then
        echo "Stuck. Retrying"
        rm $filename.jf
      fi
    done
  fi
  jellyfish dump -c $filename.jf > $filename.kmr
  rm $filename.jf
  echo "($progress/$file_count) Generate kmer count for $file [COMPLETE]"
  progress=$(($progress+1))
done

echo "Generate kmer count [COMPLETE]"
