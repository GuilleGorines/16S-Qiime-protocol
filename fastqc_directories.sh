#!/bin/bash

DIRECTORIES="$1"

for directory in $DIRECTORIES;
do
    cd $directory
    fastqc -o fastqc_result *.gz

done
