#!/bin/bash
# Regression checks for training, classification and savefiles.
# Usage: tests/run_tests.sh <NB binary>  (set LEGACY_BIN to an older build to also test format-1 savefiles)
set -u
BIN=$1; T=$(mktemp -d); cd $T
mkseq(){ python3 -c "import random,sys;random.seed($2);print('>s');print(''.join(random.choice('ACGT') for _ in range($1)))"; }
mkdata(){ # dir nclass ngen len k
  for c in $(seq 1 $2); do mkdir -p $1/C$c; for g in $(seq 1 $3); do f=$1/C$c/g$g.fasta; mkseq $4 $((c*100+g)) > $f
    jellyfish count -m $5 -s 1M -C -o $f.jf $f && jellyfish dump -c $f.jf > ${f%.fasta}.kmr && rm $f.jf; done; done; }
# savefile body tokens, skipping the "NBSAVE v k canonical" header if present
body(){ python3 -c "
t=open('$1').read().split(); print(' '.join(t[4:] if t[0]=='NBSAVE' else t))"; }
ngen(){ python3 -c "
import struct,math
v=int('$(body $1)'.split()[0]); print(round(math.exp(struct.unpack('d',struct.pack('Q',v))[0])-1))"; }
sumfreq(){ python3 -c "
import struct,math
v=int('$(body $1)'.split()[1]); print(round(math.exp(struct.unpack('d',struct.pack('Q',v))[0])))"; }
pass=0; fail=0
check(){ if [ "$2" == "$3" ]; then echo "PASS $1"; pass=$((pass+1)); else echo "FAIL $1: got '$2' want '$3'"; fail=$((fail+1)); fi; }

# 1: -n batching keeps every genome
mkdata d1 2 4 2000 5
$BIN train d1 -s s1 -k 5 -n 2 >/dev/null
check "nbatch C1 genomes" $(ngen s1/C1-save.dat) 4
check "nbatch C2 genomes" $(ngen s1/C2-save.dat) 4

# 2: -m flush on last genome of a class keeps it (~4 KB kmr files, 10 KB cap)
mkdata d2 2 3 3000 5
ls -la d2/C1/g1.kmr | awk '{print "   kmr bytes:", $5}'
$BIN train d2 -s s2 -k 5 -m 10 >/dev/null
check "memlimit C1 genomes" $(ngen s2/C1-save.dat) 3
check "memlimit C2 genomes" $(ngen s2/C2-save.dat) 3

# 3: no class trained twice in one batch, multithreaded
$BIN train d2 -s s3 -k 5 -m 10 -t 4 > l3.txt
check "no double-load" $(grep -c "double-load" l3.txt) 0
check "mt memlimit C1 genomes" $(ngen s3/C1-save.dat) 3

# 4: smoothing denominator is N + V (canonical V for k=3 is 32)
mkdir -p d4/A; printf ">s\nACGTTGCAAGGCTTAACGGATCCA\n" > d4/A/g.fasta
jellyfish count -m 3 -s 1M -C -o d4/A/g.jf d4/A/g.fasta && jellyfish dump -c d4/A/g.jf > d4/A/g.kmr
N=$(awk '{s+=$2}END{print s}' d4/A/g.kmr)
$BIN train d4 -s s4 -k 3 >/dev/null
check "sumfreq = N + V" $(sumfreq s4/A-save.dat) $((N+32))

# 5: classify with a cap smaller than one read doesn't crash with -r
$BIN classify d2 -s s2 -k 5 -m 2 -r out.csv > l5.txt; rc=$?
check "tiny memlimit + resf exit code" $rc 0
check "tiny memlimit classified all" "$(grep -c 'Genome with class' l5.txt)" 6

# 6: wrong k is rejected rather than silently mis-encoded
$BIN classify d2 -s s2 -k 6 > /dev/null 2>l6.txt; rc=$?
check "k mismatch rejected" $rc 1

# 7: k > 16 doesn't collide (two k=20 kmers equal in their last 16 bases)
mkdir -p d7/A d7/B
printf ">a\nAAAAACGTACGTACGTACGTT\n" > d7/A/g.fasta; printf ">b\nCCCCACGTACGTACGTACGTT\n" > d7/B/g.fasta
printf "AAAAACGTACGTACGTACGT 50\n" > d7/A/g.kmr; printf "CCCCACGTACGTACGTACGT 5\nGGGGGGGGGGGGGGGGGGGG 45\n" > d7/B/g.kmr
mkdir -p q7/B; printf "CCCCACGTACGTACGTACGT 1\n" > q7/B/r.kmr
$BIN train d7 -s s7 -k 20 >/dev/null
check "k=20 no collision" "$($BIN classify q7 -s s7 -k 20 | grep -o 'predicted [A-Z]*')" "predicted B"

if [ "${NEW_CHECKS:-1}" == "1" ]; then
# 8: savefile carries a header and no temp files are left behind
check "savefile header" "$(head -c 13 s2/C1-save.dat)" "NBSAVE 2 5 1"
check "no leftover tmp" "$(ls s2 | grep -c tmp)" 0

# 9: truncated savefile is an error, not a silently smaller model
cp -r s2 s9; f=s9/C1-save.dat; head -c $(( $(stat -c %s $f) / 2 )) $f > x && mv x $f
$BIN classify d2 -s s9 -k 5 > /dev/null 2> e9.txt; rc=$?
check "truncated savefile rejected" "$rc:$(grep -c truncated e9.txt)" "1:1"

# 10: savedir trained with another k / counting mode is rejected up front
$BIN classify d2 -s s2 -k 6 > /dev/null 2> e10.txt; rc=$?
check "savefile k mismatch rejected" "$rc:$(grep -c 'trained with k=5' e10.txt)" "1:1"
$BIN classify d2 -s s2 -k 5 -c 0 > /dev/null 2> e10b.txt; rc=$?
check "savefile canonical mismatch rejected" "$rc:$(grep -c 'canonical' e10b.txt)" "1:1"

# 11: classify against a missing or empty savedir fails cleanly
$BIN classify d2 -s nope -k 5 > /dev/null 2>&1; rc=$?
check "missing savedir: exit 1, not created" "$rc:$([ -e nope ] && echo created || echo absent)" "1:absent"
mkdir empty; $BIN classify d2 -s empty -k 5 > /dev/null 2>&1; rc=$?
check "empty savedir: exit 1" $rc 1

# 12: re-running training on the same genomes doesn't double count
cp -r s2 s12; $BIN train d2 -s s12 -k 5 -m 10 > l12.txt
check "retrain skips trained genomes" "$(ngen s12/C1-save.dat):$(grep -c Skipping l12.txt)" "3:6"
check "retrain leaves model unchanged" "$(cmp -s <(body s2/C1-save.dat | tr ' ' '\n' | sort) <(body s12/C1-save.dat | tr ' ' '\n' | sort) && echo same)" same

# 13: legacy (format 1) savefiles still load, with a warning
if [ -n "${LEGACY_BIN:-}" ]; then
  $LEGACY_BIN train d2 -s s13 -k 5 > /dev/null
  $BIN classify d2 -s s13 -k 5 > l13.txt 2> e13.txt; rc=$?
  check "legacy savefile loads" "$rc:$(grep -c 'Legacy savefile' e13.txt):$(grep -c 'Genome with class' l13.txt)" "0:1:6"
fi
fi

echo "passed=$pass failed=$fail"; cd /; rm -rf $T
