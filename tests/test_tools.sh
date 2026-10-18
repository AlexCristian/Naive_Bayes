#!/bin/bash
# Checks for the tools/ scripts. Usage: tests/test_tools.sh <repo root>
R=$1; T=$(mktemp -d); cd $T; pass=0; fail=0
check(){ if [ "$2" == "$3" ]; then echo "PASS $1"; pass=$((pass+1)); else echo "FAIL $1: got '$2' want '$3'"; fail=$((fail+1)); fi; }

# download(): retries failures and gives up after 5 attempts
sed -n '/^function clearLastLine/,/^}/p;/^function download/,/^}/p' $R/tools/genome_download.bash > dl.sh
cat > stub.sh <<'S'
tput(){ :; }; sleep(){ :; }
n=0; wget(){ n=$((n+1)); [ $n -ge 3 ]; }   # fails twice, then succeeds
S
check "download retries until success" "$(bash -c 'source stub.sh; source dl.sh; download u d >/dev/null; echo $n')" 3
bash -c 'source stub.sh; wget(){ echo x >> attempts; return 1; }; source dl.sh; download u d' >/dev/null 2>&1
check "download gives up after 5 attempts" "$(wc -l < attempts)" 5
check "download failure exits nonzero" "$(bash -c 'source stub.sh; wget(){ return 1; }; source dl.sh; download u d' >/dev/null 2>&1; echo $?)" 1

# nodes.dmp parsing keeps multi-word ranks intact
printf '9\t|\t32199\t|\tspecies group\t|\tX\t|\n562\t|\t561\t|\tspecies\t|\tY\t|\n' > nodes.dmp
check "nodes.dmp ranks" "$(awk -F '\t[|]\t' '{print $3}' nodes.dmp | tr '\n' ',')" "species group,species,"

# read splitter: wrapped reads stay whole, exactly N files
printf '>r1\nAAAAA\nCCCCC\n>r2\nGGGGG\nTTTTT\n>r3\nACGTA\n' > g.fasta
awk -v filename="g" -v reps=2 "$(grep -o "'BEGIN{OFS=\"\"; ORS=\"\"}[^']*'" $R/tools/generateSeqErrGens.bash | tr -d "'")" g.fasta
check "splitter keeps wrapped read whole" "$(cat g.2.fasta | tail -1)" "ggggg"ttttt
check "splitter writes reps files" "$(ls g.*.fasta | wc -l)" 2
check "grinder samples both strands" "$(grep -c -- '^ *grinder .*-un 0 ' $R/tools/generateSeqErrGens.bash)" 2

# fold split: seeded, balanced, no root class
mkdir -p .tree_data pool
printf '1\n2\n10\n11\n20\n' > .tree_data/child.txt; printf '1\n1\n2\n2\n1\n' > .tree_data/parent.txt
printf 'no rank\ngenus\nspecies\nspecies\nphylum\n' > .tree_data/class.txt
: > .tree_data/ids.txt; : > .tree_data/filenames.txt
for i in $(seq 1 10); do t=$([ $i -le 9 ] && echo $((10 + i % 2)) || echo 20); echo $t >> .tree_data/ids.txt; echo g$i.kmr >> .tree_data/filenames.txt; touch pool/g$i.kmr pool/g$i.fasta; done
cp -r pool pool_bak
python3 $R/tools/generate_dir_structure.py species pool 3 > /dev/null 2> split.err
sizes=$(for f in fold*; do find $f -name '*.kmr' | wc -l; done | sort | tr '\n' ' ')
check "folds balanced (9 genomes, 3 folds)" "$sizes" "3 3 3 "
check "genome without species skipped" "$(ls -d fold*/1 2>/dev/null | wc -l):$(grep -c 'Skipped 1' split.err)" "0:1"
a=$(cd fold1 && find . -name '*.kmr' | sort | tr '\n' ' '); rm -rf fold*; rm -rf pool; cp -r pool_bak pool
python3 $R/tools/generate_dir_structure.py species pool 3 > /dev/null 2>&1
check "fold split reproducible" "$(cd fold1 && find . -name '*.kmr' | sort | tr '\n' ' ')" "$a"

# copyFractionOfSeqs: 10 copies all, 0 copies none
mkdir -p src/c; for i in $(seq 1 20); do touch src/c/s$i.fasta; done
bash $R/tools/copyFractionOfSeqs.bash src all 10 > /dev/null; bash $R/tools/copyFractionOfSeqs.bash src none 0 > /dev/null
check "copy fraction 10/10 and 0/10" "$(ls all/c | wc -l):$(ls none 2>/dev/null | wc -l)" "20:0"

# accuracy_analysis: parses the current log format
mkdir -p results save_1
printf 'Genome with class 5, predicted 5, filename: "a/5/r1.kmr"\nGenome with class 5, predicted 7, filename: "a/5/r2.kmr"\nGenome with class 9, predicted 5, filename: "a/9/r.kmr"\n' > results/2010_j1.txt
touch save_1/5-save.dat save_1/7-save.dat
out=$(python3 $R/tools/accuracy_analysis.py --folds 1 --savedir 'save_{fold}' | tr '\n' ' ')
check "accuracy_analysis new format" "$out" "overall_accuracy: 0.3333333333333333 labeled__percent: 0.6666666666666666 labeled_accuracy: 0.5 "

echo "passed=$pass failed=$fail"; cd /; rm -rf $T
