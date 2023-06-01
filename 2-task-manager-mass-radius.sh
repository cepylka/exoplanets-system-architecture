#!/bin/bash

[ -d ./data/star_pickles ] || mkdir -p ./data/star_pickles

cat systems.txt | while read star || [[ -n $star ]];
do
   python ./2-mass-radius-only-from-NASA.py "$star"
   if [ $? -eq 0 ]; then
        echo -e "Done\n"
    else
        echo "Processing of $star failed" >&2
        exit 1
    fi
done
