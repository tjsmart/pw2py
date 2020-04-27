#!/usr/bin/env bash

for f in demo* ; do 
    echo "diff $f expected_results/$f"
    diff $f expected_results/$f
    echo
done
