#!/bin/bash
# Regression checks for training, classification and savefiles.
# Usage: tests/run_tests.sh <NB binary>  (set LEGACY_BIN to an older build to also test format-1 savefiles)
set -u
BIN=$1; T=$(mktemp -d); cd $T
mkseq(){ python3 -c "import random,sys;random.seed($2);print('>s');print(''.join(random.choice('ACGT') for _ in range($1)))"; }
mkdata(){ # dir nclass ngen len k
  for c in $(seq 1 $2); do mkdir -p $1/C$c; for g in $(seq 1 $3); do f=$1/C$c/g$g.fasta; mkseq $4 $((c*100+g)) > $f
    jellyfish count -m $5 -s 1M -C -o $f.jf $f && jellyfish dump -c $f.jf > ${f%.fasta}.kmr && rm $f.jf; done; done; }
# Savefile reader for formats 1-3 (helpers below print fields from it).
SAVEPY='
import struct,sys
def v4body(d,o):
    import struct
    ng,sf,n,nd=struct.unpack_from("<QQQQ",d,o); o+=32
    dic=struct.unpack_from("<%dQ"%nd,d,o); o+=8*nd
    w=d[o]; o+=1; (kb,)=struct.unpack_from("<Q",d,o); o+=8
    ks=[]; prev=0; p=o; e=o+kb
    while p<e:
        g=0; sh=0
        while True:
            b=d[p]; p+=1; g|=(b&0x7f)<<sh
            if not b&0x80: break
            sh+=7
        prev+=g; ks.append(prev)
    o=e; vs=[dic[int.from_bytes(d[o+i*w:o+(i+1)*w],"little")] for i in range(n)]; o+=n*w
    (m,)=struct.unpack_from("<Q",d,o); o+=8
    return ng,sf,ks,vs,sorted(d[o:].decode().split("\n")[:m])
def read(p):
    d=open(p,"rb").read()
    if not d.startswith(b"NBSAVE"):                       # format 1 (text, no header)
        t=d.split(); n=int(t[2]); kv=sorted((int(t[3+2*i]),int(t[4+2*i])) for i in range(n))
        return int(t[0]),int(t[1]),kv,[]
    nl=d.index(b"\n"); ver=int(d[:nl].split()[1])
    if ver==2:
        t=d[nl+1:].split(); n=int(t[2]); kv=sorted((int(t[3+2*i]),int(t[4+2*i])) for i in range(n))
        return int(t[0]),int(t[1]),kv,sorted(x.decode() for x in t[3+2*n+2:])
    if ver==4:
        ng,sf,ks,vs,names=v4body(d,nl+1); return ng,sf,list(zip(ks,vs)),names
    o=nl+1; ng,sf,n=struct.unpack_from("<QQQ",d,o); o+=24
    ks=struct.unpack_from("<%dQ"%n,d,o); o+=8*n; vs=struct.unpack_from("<%dQ"%n,d,o); o+=8*n
    (m,)=struct.unpack_from("<Q",d,o); o+=8
    return ng,sf,list(zip(ks,vs)),sorted(d[o:].decode().split("\n")[:m])
'
lg(){ python3 -c "$SAVEPY
import math
v=read('$1')[$2]; print(round(math.exp(struct.unpack('d',struct.pack('Q',v))[0])${3:-}))"; }
ngen(){ lg $1 0 -1; }
sumfreq(){ lg $1 1; }
model(){ python3 -c "$SAVEPY
print(read('$1'))"; }
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
check "savefile header" "$(head -c 13 s2/C1-save.dat)" "NBSAVE 4 5 1"
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
check "retrain leaves model unchanged" "$([ "$(model s2/C1-save.dat)" == "$(model s12/C1-save.dat)" ] && echo same)" same

# 14: --no-prior drops the prior term. Small (1 genome) has a slightly higher
# likelihood for the query; Big (2 genomes) only wins through its prior.
mkdir -p d14/Big d14/Small q14/Small
printf "AAAAA 1\n" > d14/Big/g1.kmr; printf "CCCCC 1\n" > d14/Big/g2.kmr; printf "AAAAA 1\n" > d14/Small/g1.kmr
printf "AAAAA 1\n" > q14/Small/r.kmr
$BIN train d14 -s s14 -k 5 > /dev/null
check "prior favors bigger class" "$($BIN classify q14 -s s14 -k 5 | grep -o 'predicted [A-Za-z]*')" "predicted Big"
check "--no-prior removes it" "$($BIN classify q14 -s s14 -k 5 --no-prior | grep -o 'predicted [A-Za-z]*')" "predicted Small"

# 15: counting k-mers straight from FASTA (gzipped or not) gives the same
# model as jellyfish count -C: several records, N runs, lowercase, short records
mkdir -p d15j/A d15f/A d15g/A
python3 -c "
import random; r=random.Random(5)
seq=lambda n: ''.join(r.choice('ACGT') for _ in range(n))
s='>r1 x\n'+seq(500)+'NNNN'+seq(300).lower()+'\n'+seq(80)+'\n>r2\n'+seq(10)+'\n>r3\n'+seq(700)+'N'+seq(40)+'\n'
open('d15j/A/g.fasta','w').write(s); open('d15f/A/g.fna','w').write(s)"
gzip -c d15f/A/g.fna > d15g/A/g.fna.gz
jellyfish count -m 11 -s 1M -C -o d15j/A/g.jf d15j/A/g.fasta && jellyfish dump -c d15j/A/g.jf > d15j/A/g.kmr
$BIN train d15j -s s15j -k 11 > /dev/null; $BIN train d15f -s s15f -k 11 -e .fna > /dev/null; $BIN train d15g -s s15g -k 11 -e .gz > /dev/null
strip(){ python3 -c "$SAVEPY
x=read('$1'); print(x[:3])"; }
check "FASTA counts = jellyfish -C" "$([ "$(strip s15j/A-save.dat)" == "$(strip s15f/A-save.dat)" ] && echo same)" same
check "gzipped FASTA counts = jellyfish -C" "$([ "$(strip s15j/A-save.dat)" == "$(strip s15g/A-save.dat)" ] && echo same)" same

# 13: legacy (format 1) savefiles still load, with a warning
if [ -n "${LEGACY_BIN:-}" ]; then
  $LEGACY_BIN train d2 -s s13 -k 5 > /dev/null
  $BIN classify d2 -s s13 -k 5 > l13.txt 2> e13.txt; rc=$?
  check "legacy savefile loads" "$rc:$(grep -c 'Legacy savefile' e13.txt):$(grep -c 'Genome with class' l13.txt)" "0:1:6"
fi
fi

echo "passed=$pass failed=$fail"; cd /; rm -rf $T
