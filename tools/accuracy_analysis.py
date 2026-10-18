#! /usr/bin/env python3
"""Per-fold accuracy from NB classify logs.

For each fold, reports accuracy over all reads, accuracy over reads whose true
class was trained on ("labeled"), and the share of such reads.
"""
import argparse
import glob
import os
import re

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--year', default='2010')
parser.add_argument('--folds', type=int, default=5)
parser.add_argument('--results', default='results/{year}_j{fold}.txt',
                    help='classify log per fold')
parser.add_argument('--savedir',
                    default='/lustre/scratch/zz374/NBC/process/temporal_data/{year}/save_{fold}',
                    help='savedir the fold was classified against')
args = parser.parse_args()

# Matches both the old line format and the one with ", filename: ..." appended.
RESULT = re.compile(r'^Genome with class (\S+), predicted ([^,\s]+)')

pref = []
for fold in range(1, args.folds + 1):
    pairs = []
    with open(args.results.format(year=args.year, fold=fold)) as f:
        for line in f:
            m = RESULT.match(line)
            if m:
                pairs.append(m.groups())

    savedir = args.savedir.format(year=args.year, fold=fold)
    # Savefiles are named <class>-save.dat; class IDs may contain '-'.
    labeled = {os.path.basename(p).rsplit('-', 1)[0]
               for p in glob.glob(os.path.join(savedir, '*.dat'))}

    known = [(t, p) for t, p in pairs if t in labeled]
    overall = sum(t == p for t, p in pairs) / len(pairs)
    labeled_acc = sum(t == p for t, p in known) / len(known) if known else float('nan')
    pref.append((overall, labeled_acc, len(known) / len(pairs)))

print('overall_accuracy:')
for overall, _, _ in pref:
    print(overall)
print('labeled__percent:')
for _, _, share in pref:
    print(share)
print('labeled_accuracy:')
for _, labeled_acc, _ in pref:
    print(labeled_acc)
