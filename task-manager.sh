#!/bin/bash

inputFolder=""
pathToScript=""
#outputFolder=""

#while getopts ":i:s:o:" opt; do
while getopts ":i:s:" opt; do
  case $opt in
    i) inputFolder="$OPTARG"
    ;;
    s) pathToScript="$OPTARG"
    ;;
    #o) outputFolder="$OPTARG"
    #;;
    \?)
        echo "[ERROR] Unknown argument -$OPTARG" >&2
        exit 1
    ;;
  esac
done

if [[
    "$inputFolder" == ""
    || "$pathToScript" == ""
    #|| "$outputFolder" == ""
]]; then
    echo "[ERROR] You need to provide input folder (-i) and path to script (-s)" >&2
    exit 2
fi

if [ ! -d $inputFolder ]
then
    echo "[ERROR] Provided input directory does not exist"
    exit 3
fi

if [ ! -f $pathToScript ]
then
    echo "[ERROR] Provided script does not exist"
    exit 4
fi

for p in $inputFolder/*.pkl
do
    python "$pathToScript" "$p"
    if [ $? -eq 0 ]
    then
        echo -e "Done\n"
    else
        echo "Processing of $p failed" >&2
        exit 1
    fi
done
