#! /bin/bash

function clearLastLine {
  tput cuu 1 && tput el
}

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]
  then
    echo "usage: ./reorganize.bash [taxonomic rank] [# folds] [source folder to be reorganized]
    kingdom
    class
    subclass
    superorder
    order
    suborder
    family
    subfamily
    genus
    species
    "
    exit
  else
    DEPTH=$1
fi

TARGET_DIR=$3
POOL_FOLDER="$TARGET_DIR/pool"
N_FOLDS=$2


if [ -d $POOL_FOLDER ] || [ -f $POOL_FOLDER ]
then
  rm -r $POOL_FOLDER
fi
mkdir $POOL_FOLDER

echo "Move all genomes to pool folder [IN PROGRESS]"

for dir in $(ls $TARGET_DIR)
do
  dirpath=$TARGET_DIR/$dir
  if [ $dirpath != $POOL_FOLDER ]
  then
    echo "Looking at $dirpath"
    find $dirpath -iname '*.kmr' -exec mv {} $POOL_FOLDER \;
    find $dirpath -iname '*.fasta' -exec mv {} $POOL_FOLDER \;
    rm -r $dirpath
  fi
done
#clearLastLine
echo "Move all genomes to pool folder [COMPLETE]"


awk 'BEGIN{FS=OFS="/";filesuffix="genomic.kmr"}{id=$1;sub(" ftp:", "", id);asm=$10;file=asm"_"filesuffix;print id" "file}' $TARGET_DIR/.tree_data/ftpdirpaths > $TARGET_DIR/.tree_data/id_filename.txt
awk '{print $1}' $TARGET_DIR/.tree_data/id_filename.txt > $TARGET_DIR/.tree_data/ids.txt
awk '{print $2}' $TARGET_DIR/.tree_data/id_filename.txt > $TARGET_DIR/.tree_data/filenames.txt

chmod +x generate_dir_structure.py
./generate_dir_structure.py $DEPTH $POOL_FOLDER $N_FOLDS
#rm -r $POOL_FOLDER

echo "Done."
