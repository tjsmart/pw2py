#!/usr/bin/env bash

delete="New_Si.sample.in  New_Si.sample.out input_tmp.in temp CRASH"
for dlt in delete; do
    rm -rf $dlt 2> /dev/null
done
echo "Clean"
