#! /bin/bash

if [ -z "$1" ]
  then
    echo "usage: ./genome_download.bash [taxonomic rank]
    domain
    kingdom
    division
    subdivision
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
    depth=$1
fi

function clearLastLine {
  tput cuu 1 && tput el
}

function download {
  if [ -z "$3" ]
  then
    down_text="Download"
  else
    down_text="($3) Download"
  fi

  echo "$down_text $1 - [IN PROGRESS]"
  wget $1 --directory-prefix $2 -q && retval=$?
  while [ $retval -ne "0" ]
  do
    clearLastLine
    echo "$down_text $1 - [FAILED]"
    echo "Retrying."
    echo "$down_text $1 - [IN PROGRESS]"
    wget $1 --directory-prefix $2 -q && retval=$?
  done

  clearLastLine
  echo "$down_text $1 - [COMPLETE]"
}

IF_DOWNLOAD_TREEDATA=$2

if [ $IF_DOWNLOAD_TREEDATA == "true" ]
then

    if [ -d .tree_data ]
    then
      rm -r .tree_data
    fi
    mkdir .tree_data # Create a throwaway dir for auxiliary files

    echo "Download genome listings [IN PROGRESS]"

    # These files change as more genomes are added
    download ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt .tree_data # Correlates IDs with FTP download links
    download ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz .tree_data # Correlates node IDs with parent IDs
    tar -zxf .tree_data/taxdump.tar.gz -C .tree_data nodes.dmp # Extract only what we need
else
    echo "Constructing from existing tree data [IN PROGRESS]"
fi
# Extract the parents and children from the node file, and split the data
# over two files.
awk '{print $1}' .tree_data/nodes.dmp > .tree_data/child.txt
awk '{print $3}' .tree_data/nodes.dmp > .tree_data/parent.txt
awk '{print $5}' .tree_data/nodes.dmp > .tree_data/class.txt

# Begin processing data into FTP links
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $6, $20}' .tree_data/assembly_summary.txt > .tree_data/ftpdirpaths # Preserves ID and starts printing FTP link

# Output in the format "ID ftpLink"
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' .tree_data/ftpdirpaths > .tree_data/ftpfilepaths
awk '{print $1}' .tree_data/ftpfilepaths > .tree_data/ids.txt
awk '{print $2}' .tree_data/ftpfilepaths > .tree_data/ftp_links.txt

# Finished building tree, start building directory structure
if [ -d genomes ]
then
  rm -r genomes
fi
mkdir genomes

echo "Download genome listings [COMPLETE]"
echo "Constructing specified tree [IN PROGRESS]"

chmod +x generate_dir_structure.py
./generate_dir_structure.py $depth genomes

clearLastLine
echo "Constructing specified tree [COMPLETE]"
echo "Download genomes [IN PROGRESS]"

file_count=$(cat .tree_data/ftpdirpaths | wc -l)
progress=1

# Proceed to download the files
for dir in $(ls -d genomes/*/)
do
  download_list=$dir"to_download.txt"
  while read url
  do
    download $url $dir $progress/$file_count
    progress=$(($progress+1))
  done < $download_list
  rm $download_list
done

clearLastLine
echo "Download genomes [COMPLETE]"
echo "Decompress genomes [IN PROGRESS]"

# Decompress the files
progress=1
for archive in $(ls genomes/*/*.gz)
do
  echo "($progress/$file_count) Decompress $archive [IN PROGRESS]"
  gzip -dq $archive
  clearLastLine
  echo "($progress/$file_count) Decompress $archive [COMPLETE]"
  progress=$(($progress+1))
done

echo "Decompress genomes [COMPLETE]"

clearLastLine
echo "Download genomes [COMPLETE]"

#./gbfftofasta.bash genomes

echo "Download successful."
