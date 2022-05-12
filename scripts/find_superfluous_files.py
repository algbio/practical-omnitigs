#!/usr/bin/env python3

import os, sys

with open("/wrk-vakka/users/sebschmi/practical-omnitigs/all_files.txt", "r") as all_files_txt:
    all_files = set([f.strip() for f in all_files_txt.readlines()])

with open("/proj/sebschmi/git/practical-omnitigs/required_files.txt", "r") as required_files_txt:
    required_files = set([f.strip() for f in required_files_txt.readlines()])

if "" in all_files:
    all_files.remove("")

if "" in required_files:
    required_files.remove("")

valid_prefixes = set()
for required_file in required_files:
    valid_prefixes.add(os.path.dirname(required_file))

all_prefixes = set()
for f in all_files:
    all_prefixes.add(os.path.dirname(f))

#for valid_prefix in sorted(valid_prefixes):
#    print(valid_prefix)

superfluous_files = set()

for f in all_files:
    superfluous = True
    for valid_prefix in valid_prefixes:
        if f.startswith(valid_prefix):
            superfluous = False
            break
        else:
            assert not valid_prefix.startswith(f)

    if superfluous and f not in superfluous_files:
        #print(f)
        superfluous_files.add(f)

for superfluous_file in sorted(superfluous_files):
    print(superfluous_file)
