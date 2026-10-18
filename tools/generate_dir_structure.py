#! /usr/bin/env python3
import os
from os.path import splitext, dirname
import sys
from enum import Enum
import random

__author__ = "Alexandru Cristian"

'''
    Internal script - should not be called outside of reorganize.bash
                        ~Programmer's reference~
    Argument list:
        ./generate_dir_structure.py (1) [taxonomic rank, string]
                                    (2) [pool folder (source dir where all files are dumped)]
                                    (3) [number of folds]
'''

'''
    DEPRECATED:
    Internal script - should not be called outside of reorganize.bash
                        ~Programmer's reference~
    Argument list:
        ./generate_dir_structure.py (1) [taxonomic rank, string]
                                    (2) [pool folder (source dir where all files are dumped)]
                                    (3) [target directory for training set]
                                    (4) [target directory for testing set]
                                    (5) [distribution split factor (0-10)]
'''

class GenMode(Enum):
    GEN_FTP_LINKS = 1
    REGEN_FROM_POOL = 2

mode = GenMode.GEN_FTP_LINKS
if len(sys.argv) > 3:
    mode = GenMode.REGEN_FROM_POOL

ORIGIN = "."
POOL_FOLDER = "."
FOLD_BASENAME = ""
if mode is GenMode.REGEN_FROM_POOL:
    ORIGIN = dirname(sys.argv[2])
    if ORIGIN == "": # Fix Python's dirname inconsistency with bash's response
        ORIGIN = "."
    POOL_FOLDER = sys.argv[2]
    FOLD_BASENAME = ORIGIN + "/fold"
    N_FOLDS = int(sys.argv[3])
    # Seeded, so the same inputs always give the same split (override with a
    # 4th argument).
    SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 0

tree = dict()
with open(ORIGIN + "/.tree_data/child.txt") as children, open(ORIGIN + "/.tree_data/parent.txt") as parents, open(ORIGIN + "/.tree_data/class.txt") as classes:
    for child, parent, cls in zip(children, parents, classes):
        tree[child.strip()] = (parent.strip(), cls.strip())

depth = sys.argv[1]
if mode == GenMode.GEN_FTP_LINKS:
    urls_file = ORIGIN + "/.tree_data/ftp_links.txt"
else:
    urls_file = ORIGIN + "/.tree_data/filenames.txt"

with open(ORIGIN + "/.tree_data/ids.txt") as node_ids, open(urls_file) as urls:
    entries = []
    skipped = 0
    for node_id, url in zip(node_ids, urls):
        node_id = node_id.strip()
        url = url.strip()
        while tree[node_id][1] != depth and node_id != "1":
            node_id = tree[node_id][0]
        if node_id == "1":
            # No ancestor at this rank: don't lump these genomes together
            # into one bogus class named after the taxonomy root.
            skipped += 1
            continue
        entries.append((node_id, url))
    if skipped:
        print("Skipped %d genomes with no %s-level ancestor" % (skipped, depth),
              file=sys.stderr)

if mode is GenMode.REGEN_FROM_POOL:
    # Shuffle, then deal genomes round-robin so fold sizes differ by at most one.
    random.Random(SEED).shuffle(entries)

for n, (node_id, url) in enumerate(entries):
    if mode is GenMode.REGEN_FROM_POOL:
        path = FOLD_BASENAME + str(n % N_FOLDS + 1) + "/" + node_id
    elif mode is GenMode.GEN_FTP_LINKS:
        path = ORIGIN + "/" + sys.argv[2] + "/" + node_id
    print(path)
    if not os.path.isdir(path):
        os.makedirs(path)
    if mode is GenMode.GEN_FTP_LINKS:
        with open(path + '/to_download.txt', 'a', encoding='utf-8') as file:
            file.write(url + "\n")
    else:
        os.rename(POOL_FOLDER + "/" + url, path + "/" + url)
        seq_url = splitext(url)[0] + ".fasta"
        os.rename(POOL_FOLDER + "/" + seq_url, path + "/" + seq_url)
